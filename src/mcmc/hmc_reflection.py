#!/usr/bin/env python
# coding: utf-8

#################################################
# This module provides a basic HMC sampler,
#   which will be used as a template to develope
#   more samplers.
#
# This code is written by:
#   Ahmed Attia (attia@anl.gov)
#   2023
# To use all or parts of the code,
#   you must get concent from the Author
#################################################

## Import
# Essentials for Class design
from abc import (
    ABC as _ABC,
    abstractmethod as _abstractmethod,
)

# Math & Science
import numbers
import numpy as np

try:
    from scipy import sparse
    from scipy.sparse import csc_matrix
    from scipy.sparse import linalg
except ImportError:
    sparse = None

# Other imports
import time
import re

import multiprocess
from concurrent.futures import ProcessPoolExecutor

## Module-level variables
_DEBUG = False

print("Importing MCMC Module: mcmc_sampling_reflection.py")


class Sampler(_ABC):
    def __init__(self, configs={}):
        """
        Base class for Samplers (algorithms to generate samples from a predefined distribution).
        An example is an inverse CDF sampler, MCMC sampler, etc.

        :param dict configs: a key-value dictionary that holds the configurations of the sampler

        :raises:
            - :py:class:`TypeError` is raised if the passed `configs` is not a dictionary
            - :py:class:`TypeError` is raised if the passed `configs` holds invalid keys/values/types

        :remarks:
            - This method is either copied into derived methods, or is called (by using super().__init__ properly)
        """
        if not isinstance(configs, dict):
            print(
                f"The passed configs object must be a dictionary; received {type(configs)}"
            )
            raise TypeError
        else:
            if self.validate_configurations(configs):
                self._CONFIGURATIONS = configs
            else:
                print(
                    f"The passed configurations dictionary includes invlid values/types"
                )
                raise TypeError

    @_abstractmethod
    def sample(self, sample_size):
        """
        This function generates samples from the created sampler.
        This method returns a list with each entry representing a sample point from the underlying distribution
        """
        ...

    @_abstractmethod
    def validate_configurations(self, configs):
        """
        A method to check the passed configuratios and make sure they
            are conformable with each other, and with current configurations (or default of not set)

        :returns:
            True/False flag indicating whether passed coinfigurations dictionary is valid or not
        """
        ...

    def update_configurations(self, **kwargs):
        """
        Take any set of keyword arguments, and lookup each in
        the configurations, and update as nessesary/possible/valid

        :raises:
            - :py:class:`TypeError` is raised if any of the passed keys in `kwargs` is invalid/unrecognized

        :remarks:
            - Generally, we don't want actual implementations in abstract classes, however, this one is provided
            as good guidance. Derived classes can rewrite it.
        """
        if self.validate_configurations(kwargs, raise_for_invlid=True):
            # Update underlying settings
            self._CONFIGURATIONS.update(kwargs)

    def _aggregate_configurations(self, configs, default_configs):
        """
        Given two dictionaries `configs`, `default_configs`, with the latter holding default values,
            aggregate/combine the key-value pairs in the two dictionaries with precedence given to the configs.
            Only keys those are exclusively in `default_configs`` are taken, and the common ones are taken from `configs`.

        :raises:
            - :py:class:`TypeError` is raised if the passed `configs` is not a dictionary

        :remarks:
            - Generally, we don't want actual implementations in abstract classes, however, this one is provided
            as good guidance. Derived classes can rewrite it.
        """
        if not isinstance(configs, dict) or not isinstance(default_configs, dict):
            print(
                f"Both of the two configs objects must be a dictionary; received {type(configs)}, {type(default_configs)}"
            )
            raise TypeError

        updated_configs = default_configs.copy()
        for key, value in configs.items():
            if isinstance(key, dict):
                updated_configs[key] = self._aggregate_configurations(
                    value, updated_configs[key]
                )
            else:
                updated_configs[key] = value
        return updated_configs

    @property
    def configurations(self):
        """
        A property to return a copy (not the dictionary itself) of the underlying configurations.
        To change configurations, call :py:meth:`update_configurations`
        """
        return self._CONFIGURATIONS.copy()


class FDGradient(object):
    def __init__(self, func, x, fd_eps=1e-5):
        """Create a gradient object"""
        self.x = np.asarray(x, dtype=float).flatten()
        self.e = np.zeros_like(self.x)
        self.func = lambda x: float(func(x))
        self.fd_eps = fd_eps

    def __call__(self, i):
        """Evaluate and return ith entry of the gradient"""
        self.e[...] = 0
        self.e[i] = self.fd_eps
        fd = float(self.func(self.x + self.e) - self.func(self.x)) / self.fd_eps
        return fd


class HMCSampler(Sampler):
    # A dictionary holding default configurations
    _DEF_CONFIGURATIONS = {
        "size": None,
        "log_density": None,
        "log_density_grad": None,
        "parallel_fd_grad": False,
        "random_seed": None,
        "burn_in": 500,
        "mix_in": 10,
        "symplectic_integrator": "verlet",
        "symplectic_integrator_stepsize": 1e-2,
        "symplectic_integrator_num_steps": 20,
        # 'mass_matrix':1,
        "constraint_test": None,
        "mass_matrix_gamma": np.ones(25),
        "mass_matrix_theta": np.ones(25),
        "SD_gamma": np.ones(25),
        "SD_theta": np.ones(25),
    }

    def __init__(self, configs=_DEF_CONFIGURATIONS):
        """
        Implementation of the HMC sampling algorithm
            (with multiple choices of the symplectic integrator).

        :param dict configs: a configurations dictionary wchich accepts the following keys:
            - 'size': (int) dimension of the target distribution to sample
            - 'log_density': (callable) log of the (unscaled) density function to be sampled
            - 'log_density_grad': (callable) the gradient of the `log_density` function passed.
                If None, an attempt will be made to used automatic differentiation
                (if available), otherwise finite differences (FD) will be utilized
            - parallel_fd_grad: if `True` evaluate the FD gradient in parallel (if FD is used,
                that is when log_density_grad is not passed)
            - random_seed: random seed used when the object is initiated to keep track of random samples
              This is useful for reproductivity.
              If `None`, random seed follows `numpy.random.seed` rules
            - 'burn_in': (int) number of sample points to discard before collecting samples
            - 'mix_in': (int) number of generated samples between accepted ones
                (to descrease autocorrelation)
            - 'symplectic_integrator': (str) name of the symplectic integrator to use;
                acceptable are:
                  + 'leapfrog', '2-stage', '3-stage',
                      where both 'leapfrog' and 'verlet' are equivalent
            'symplectic_integrator_stepsize': (positive scalar) the step size of the symplectic
                integrator
            'symplectic_integrator_num_steps': (postive integer) number of steps of size
                `symplectic_integrator_stesize` taken before returnig the next proposed point
                over the Hamiltonian trajectory
            'mass_matrix': (nonzero scalar or SPD array of size `size x size`)  mass matrix
                to be used to adjust sampling the auxilliary Gaussian momentum
            'constraint_test': a function that returns a boolean value `True` if sample point satisfy any constrints, and `False` otherwise

        :remarks:
            - Validation of the configurations dictionary is taken care of in the super class
            - References:
        """
        configs = self._aggregate_configurations(configs, self._DEF_CONFIGURATIONS)
        super().__init__(configs)

        # Define additional private parameters
        self._update_mass_matrix()

        # Define log-density and associated gradient
        self._update_log_density()

        # maintain a proper random state
        random_seed = self._CONFIGURATIONS["random_seed"]
        self._RANDOM_STATE = np.random.RandomState(random_seed).get_state()

    def __factorize_spsd_matrix(self, A, shape=None, lower=True):
        """
        Genereate Cholesky decomposition of the passed symmetric positive semidefinite matrix.

        :param A: `float`, or 2d `numpy.ndarray` or `scipy.spmatrix`

        :returns: square root of `A` if scalar, otherwise, lower triangular portion of
            the matrix Cholesky decomposition.

        :remarks:
            This function is a utility function and should be moved to a utility module as it is indenpendent from this class
        """
        if isinstance(A, numbers.Number):
            std = np.sqrt(A)

        elif isinstance(A, np.ndarray):
            n = A.shape[0]
            assert A.ndim == 2 and A.shape == (
                n,
                n,
            ), "The matrix must by square of 2dimensions"

            # Eliminate rows/columns with all zeros
            valid_inds = np.asarray(
                [
                    i
                    for i in range(n)
                    if not (np.count_nonzero(A[i, :]) == np.count_nonzero(A[:, i]) == 0)
                ]
            )
            zero_inds = np.setdiff1d(np.arange(n), valid_inds)

            if valid_inds.size == 0:
                std = np.zeros_like(A, dtype=A.dtype)
            else:
                A_reduced = A[valid_inds, :][:, valid_inds]

                # Cholesky factor (for nonzero rows/columns)
                std_reduced = np.linalg.cholesky(A_reduced)
                if not lower:
                    std_reduced = std_reduced.T

                std = std_reduced
                if zero_inds.size > 0:
                    for i in zero_inds:
                        for axis in [0, 1]:
                            std = np.insert(std, i, 0.0, axis=axis)

        elif sparse.issparse(A):
            n = A.shape[0]
            assert A.ndim == 2 and A.shape == (
                n,
                n,
            ), "The matrix must by square of 2dimensions"

            # make sure the maitrix is indexible
            if A.format in ["dia", "coo"]:
                A = A.tocsc()

            # Eliminate rows/columns with all zeros
            valid_inds = np.asarray(
                [i for i in range(n) if not (A[[i], :].nnz == A[:, [i]].nnz == 0)]
            )
            zero_inds = np.setdiff1d(np.arange(n), valid_inds)

            if valid_inds.size == 0:
                std = sparse.csc_matrix((n, n))
            else:
                A_reduced = A[valid_inds, :][:, valid_inds]

                # Sparse LU (for the reduced version)
                LU = sparse.linalg.splu(A_reduced, diag_pivot_thresh=0)

                # Check the matrix is positive semi definite:
                if np.all(LU.perm_r == np.arange(valid_inds.size)) and np.all(
                    LU.U.diagonal() > 0
                ):
                    std_reduced = LU.L.dot(sparse.diags(LU.U.diagonal() ** 0.5))
                    if not lower:
                        std_reduced = std_reduced.T
                    std_reduced = std_reduced.tocsc()

                    # Insert zeros into std.
                    if zero_inds.size == 0:
                        std = std_reduced
                    else:
                        # Expand
                        std = sparse.lil_matrix((n, n))
                        for i in range(valid_inds.size):
                            std[valid_inds[i], valid_inds] = std_reduced[i, :]
                            # std[:, valid_inds[i]] = std_reduced[i]

                        std = std.tocsc()

                else:
                    print(
                        "Failed to use efficient LU factorization for sparse matrices;"
                    )
                    print("Trying dense!")
                    std = factorize_spsd_matrix(A.toarray(), shape=shape, lower=lower)
                    std = sparse.csc_matrix(std)

        else:
            print("A must be scalar, numpy array,scipy sparse matrix ")
            print(f"Received {type(A)}")
            raise TypeError

        return std

    def _update_mass_matrix(self):
        """
        Update the momentum covariance, i.e., the mass matrix given the current mass matrix
            in the configurations dictionary, and the iverse of the mass matrix

        This method defines/updates three variables:
            `_MASS_MATRIX` and `_MASS_MATRIX_INV`, `_MASS_MATRIX_SQRT`  which should never be updated manually
        """
        size = self._CONFIGURATIONS["size"]
        # mass_matrix = self._CONFIGURATIONS['mass_matrix']
        mass_matrix_gamma = self._CONFIGURATIONS["mass_matrix_gamma"]
        mass_matrix_theta = self._CONFIGURATIONS["mass_matrix_theta"]

        #         if isinstance(mass_matrix, numbers.Number):
        # #             if mass_matrix <= 0:
        # #                 print(f"NonPositive momentum covariance (mass) value!")
        # #                 raise ValueError

        if sparse is not None:
            # mass_matrix = sparse.diags(np.block([[mass_matrix[0]*np.ones(int(size/2)),mass_matrix[1]*np.ones(int(size/2))]]),[0])
            # first_block = mass_matrix[0]*np.ones(int(size/2))
            # second_block = mass_matrix[1]*np.ones(int(size/2))

            # modify the first two elements of the first block
            # mass_matrix_theta[0] *= 100
            # mass_matrix_theta[1] *= 100
            # mass_matrix_theta[2] *= 10
            # mass_matrix_theta[18] *= 10
            # print("scales theta")
            # mass_matrix = sparse.diags(np.block([[first_block,second_block]]),[0])

            print(
                "mass matrix used in HMC:",
                np.block([[mass_matrix_theta, mass_matrix_gamma]]),
            )
            mass_matrix = sparse.diags(
                np.block([[mass_matrix_theta, mass_matrix_gamma]]), [0]
            )

        else:
            mass_matrix = np.array([[mass_matrix]])

        #         elif isinstance(mass_matrix, np.ndarray):
        #             if mass_matrix.shape != (size, size):
        #                 print(f"The mass matrix found has wrong shape"
        #                       f" > Expected: matrix/array of shape ({size}, {size})"
        #                       f" > Found matrix of shape: {mass_matrix.shape}")
        #                 raise TypeError

        #         elif sparse is not None and isinstance(mass_matrix, sparse.spmatrix):
        #             if mass_matrix.shape != (size, size):
        #                 print(f"The mass matrix has wrong shape of {mass_matrix.shape}"
        #                       f" > Expected: matrix/array of shape ({size}, {size})")
        #                 raise TypeError

        #         else:
        #             print(f"Invalid type of the mass matrix {type(mass_matrix)}!"
        #                   f"Expected, scalar, np array or sparse matrix/array")
        #             raise TypeError

        # Create the inverse of the mass matrix once.
        if sparse is not None:
            self._MASS_MATRIX = csc_matrix(mass_matrix)
            self._MASS_MATRIX_INV = sparse.linalg.inv(self._MASS_MATRIX)

        else:
            self._MASS_MATRIX = np.asarray(mass_matrix)
            self._MASS_MATRIX_INV = np.linalg.inv(self._MASS_MATRIX)

        # Mass matrix square root (lower Cholesky factor for sampling)
        self._MASS_MATRIX_SQRT = self.__factorize_spsd_matrix(self._MASS_MATRIX)

    def __parallel_func_grad(
        self,
        x,
        processes=multiprocess.cpu_count(),
    ):
        """Evaluate the"""
        func = self._LOG_DENSITY
        evaluator = FDGradient(
            func=func,
            x=x,
        )
        with multiprocess.Pool(processes) as pool:
            grad = pool.map(evaluator, range(len(x)))
            return np.asarray(grad)

    def __threaded_func_grad(
        self,
        x,
        processes=multiprocess.cpu_count(),
    ):
        """Evaluate the"""
        func = self._LOG_DENSITY
        evaluator = FDGradient(
            func=func,
            x=x,
        )
        with ProcessPoolExecutor() as executor:
            print("Evaluating Gradient with Multithreading")
            grad = np.asarray([v for v in executor.map(evaluator, range(len(x)))])
            print("Done.")
            return grad

    def __fd_grad_entry(self, x, i):
        e = np.zeros_like(x)
        e[i] = fd_eps

        grad_i = (func(x + e) - func(x)) / fd_eps
        return grad_i

    def __create_threaded_func_grad(
        self,
        processes=multiprocess.cpu_count(),
    ):
        """Evaluate the"""

        def func_grad(
            x,
            fd_eps=1e-5,
        ):
            """Function to generate gradient using finite differences"""
            x = np.asarray(x).flatten()
            xl = [x for _ in range(x.size)]

            with ProcessPoolExecutor(processes) as executor:
                print("Evaluating Gradient with Multithreading")
                grad = np.asarray(
                    [
                        v
                        for v in executor.map(
                            self.__fd_grad_entry, zip(xl, range(x.size))
                        )
                    ]
                )
                print("Done")

        return func_grad

    def __create_func_grad(
        self, func, size, approach="fd", fd_eps=1e-5, fd_central=False
    ):
        """
        Given a callable/function `func`, create a function that evaluates the gradient of this function
        :param int size: the domain size which determines the size of the returned gradient
        :param str approach: the approach to use for creating the function.

        :remarks: this method is planned to enable automatic differentiation (AD) if available on
            the current platform, otherwise finite differences 'fd' is used
        """
        if re.match(
            r"\A(f(-|_| )*d|finite(-|_| )*difference(s)*)\Z", approach, re.IGNORECASE
        ):

            def func_grad(x, fd_eps=fd_eps, fd_central=fd_central):
                """Function to generate gradient using finite differences"""
                x = np.asarray(x).flatten()
                grad = np.zeros_like(x)
                e = np.zeros_like(x)
                for i in range(e.size):
                    e[:] = 0.0
                    e[i] = fd_eps

                    if fd_central:
                        grad[i] = (func(x + e) - func(x - e)) / (2.0 * fd_eps)
                    else:
                        grad[i] = (func(x + e) - func(x)) / fd_eps
                return grad

            return func_grad

        elif re.match(
            r"\A(a(-|_| )*d|automatic(-|_| )*differentiation)\Z",
            approach,
            re.IGNORECASE,
        ):
            raise NotImplementedError(
                "TODO: A/D is not yet supported for creating function gradient"
            )

        else:
            print(f"Unrecognized gradient generation approach {approach}")
            raise ValueError

    def _update_log_density(self):
        """
        Update the function that evaluates the logarithm of the (unscaled) target density function
            and the associated gradient (if given). If the gradient is not given,
            either automatic differentiation (if installed/requested) is utilized, otherwise
            finite-differences are used

        This method defines/updates two variables:
            `_LOG_DENSITY` and `_LOG_DENSITY_GRAD` which evalute the value and the gradient of
                the log-density function of the (unscaled) target distribution
        """
        size = self._CONFIGURATIONS["size"]

        # Log-Density function
        log_density = self._CONFIGURATIONS["log_density"]
        if not callable(log_density):
            print(
                f"The 'log_density' found in the configurations is not a valid callable/function!"
            )
            raise TypeError
        try:
            test_vec = np.random.randn(size)
            log_density(test_vec)
        except:
            print(
                f"Failed to evaluate the log-density using a randomly generated vector"
            )
            raise TypeError
        # Associate to self
        self._LOG_DENSITY = log_density

        # Log-Density gradient
        log_density_grad = self._CONFIGURATIONS["log_density_grad"]
        if log_density_grad is None:
            # Serial version of the gradient evaluation using FD
            log_density_grad_serial = self.__create_func_grad(log_density, size=size)
            # Parallel versions of the gradient evaluation using FD
            log_density_grad_parallel = lambda x: self.__parallel_func_grad(
                x, processes=min(size, multiprocess.cpu_count())
            )
            log_density_grad_threaded = lambda x: self.__threaded_func_grad(
                x, processes=min(size, multiprocess.cpu_count())
            )
            # log_density_grad = self.__create_threaded_func_grad(processes=min(size, multiprocess.cpu_count()))

            # Test serial gradient
            test_vec = np.random.randn(size)

            if self._CONFIGURATIONS["parallel_fd_grad"]:
                print("Testing parallel gradient evaluation")
                # Test Serial vs. Paralle gradient timing
                parallel_failed = threaded_failed = True  # Initialization

                try:
                    start_time = time.time()
                    _ = log_density_grad_parallel(test_vec)
                    print(
                        f"Parallel gradient (with multiprocessing) took: {time.time()-start_time} seconds"
                    )
                    parallel_failed = False
                except:
                    print("Failed to use Parallel Gradient with MultiProcessing")

                # ONLY if Multiprocess version failed try multithreading
                if parallel_failed:
                    try:
                        start_time = time.time()
                        _ = log_density_grad_threaded(test_vec)
                        print(
                            f"Threaded gradient (with multithreading) took: {time.time()-start_time} seconds"
                        )
                        threaded_failed = False
                    except:
                        print("Failed to use Parallel Gradient Multithreading")

                if not (parallel_failed and threaded_failed):
                    # Test Serial Gradient timing
                    start_time = time.time()
                    _ = log_density_grad_serial(test_vec)
                    print(f"Serial gradient took: {time.time()-start_time} seconds")

                if parallel_failed and threaded_failed:
                    # revert to serial
                    print(
                        "Cannot generate gradient in parallel.\n"
                        "Neither The Prallel Nor the Threaded gradient could be executed.\n"
                        "Reverting to serial gradient evaluation"
                    )
                    self._CONFIGURATIONS["parallel_fd_grad"] = False

                    log_density_grad = log_density_grad_serial

                elif not parallel_failed:
                    log_density_grad = log_density_grad_parallel

                elif not threaded_failed:
                    log_density_grad = log_density_grad_threaded

                else:
                    print("This is not possible to show; report this bug!")
                    raise ValueError

            else:
                log_density_grad = log_density_grad_serial

        elif not callable(log_density_grad):
            print(
                f"The 'log_density_grad' found in the configurations is not a valid callable/function!"
            )
            raise TypeError

        try:
            test_vec = np.random.randn(size)
            grad = log_density_grad(test_vec)
            if grad.size != test_vec.size:
                print("The log density function returns gradient of wront shape")
                print("Expected gradient of size {test_vec.size}")
                print("Received gradient of size {grad.size}")
                raise AssertionError
        except:
            print(
                f"Failed to evaluate the log-density using a randomly generated vector"
            )
            raise TypeError
        # Associate to self
        self._LOG_DENSITY_GRAD = log_density_grad

    def validate_configurations(self, configs, raise_for_invalid=True):
        """
        A method to check the passed configuratios and make sure they
            are conformable with each other, and with current configurations once combined.
        This guarantees that any key-value pair passed in configs can be properly used

        :param bool raise_for_invalid: if `True` raise :py:class:`TypeError` for invalid configrations key

        :returns:
            True/False flag indicating whether passed coinfigurations dictionary is valid or not

        :raises: see the parameter `raise_for_invalid`
        """
        # Set initial value of the validity flag
        is_valid = True

        ## Check that all passed configurations keys are acceptable
        # Get a copy of current configurations, and check all passed arguments are valid
        try:
            current_configs = self._CONFIGURATIONS
        except AttributeError:
            current_configs = self.__class__._DEF_CONFIGURATIONS

        valid_keys = current_configs.keys()
        for key in configs.keys():
            if key not in valid_keys:
                if raise_for_invalid:
                    print(f"Invalid configurations key {key} passed!")
                    raise TypeError
                else:
                    is_valid = False
                    return is_valid

        # Dict containing all configurations (aggregated with current settings)
        aggr_configs = self._aggregate_configurations(configs, current_configs)

        ## Check the target distribution dimensionality
        size = aggr_configs["size"]
        if not (isinstance(size, numbers.Number) and int(size) == size and size > 0):
            if raise_for_invalid:
                print(
                    f"The size key must be a valid positive integer value; "
                    f"received {size} of type {type(size)}"
                )
                raise TypeError
            else:
                is_valid = False
                return is_valid

        ## Check the log_density function
        log_density = aggr_configs["log_density"]
        if not callable(log_density):
            if raise_for_invalid:
                print(f"The 'log_density' is not a valid callable/function")
                raise TypeError
            else:
                is_valid = False
                return is_valid

        # test ability to properly call the log_density function
        try:
            test_vec = np.random.randn(size)
            log_density(test_vec)
        except:
            if raise_for_invalid:
                print(
                    f"Failed to evaluate the log-density using a randomly generated vector"
                )
                raise TypeError
            else:
                is_valid = False
                return is_valid

        # Test constraint function (if provided)
        constraint_test = aggr_configs["constraint_test"]
        if callable(constraint_test):
            try:
                assert constraint_test(np.random.rand(size)) in [
                    False,
                    True,
                ], "Constraint test must report True/False!"
            except:
                print("The constraint test didn't work as expected!")
                raise
        else:
            assert (
                constraint_test is None
            ), "`constraint_test` must be either a callable of None!"

        ## Mass matrix (covariance of the momentum)
        # mass_matrix = aggr_configs['mass_matrix']

        #         if isinstance(mass_matrix, numbers.Number):
        #             if mass_matrix <= 0:
        #                 if raise_for_invalid:
        #                     print(f"NonPositive momentum covariance (mass) value!")
        #                     raise TypeError
        #                 else:
        #                     is_valid = False
        #                     return is_valid

        #         elif isinstance(mass_matrix, np.ndarray):
        #             if mass_matrix.shape != (size, size):
        #                 if raise_for_invalid:
        #                     print(f"The mass matrix found has wrong shape"
        #                           f" > Expected: matrix/array of shape ({size}, {size})"
        #                           f" > Found matrix of shape: {mass_matrix.shape}")
        #                     raise TypeError
        #                 else:
        #                     is_valid = False
        #                     return is_valid

        #         elif sparse is not None and isinstance(mass_matrix, sparse.spmatrix):
        #             if mass_matrix.shape != (size, size):
        #                 print(f"The mass matrix has wrong shape of {mass_matrix.shape}"
        #                       f" > Expected: matrix/array of shape ({size}, {size})")
        #                 raise TypeError

        #         else:
        #             print(f"Invalid type of the mass matrix {type(mass_matrix)}!"
        #                   f"Expected, scalar, np array or sparse matrix/array")
        #             raise TypeError

        # TODO: Proceed here
        ## Test the Symplectic symplectic_integrator parameters (symplectic_integrator name, step size, number of steps)
        pass  # TODO:

        # Return the validity falg if this point is reached (all clear)
        return is_valid

    def sample(
        self,
        sample_size=1,
        verbose=False,
        initial_state=None,
    ):
        """
        Generate and return a sample of size `sample_size`.
        This method returns a list with each entry representing a sample point from the underlying distribution
        """
        hmc_results = self.start_MCMC_sampling(
            sample_size=sample_size,
            verbose=verbose,
            initial_state=initial_state,
        )
        return hmc_results["collected_ensemble"]

    def map_estimate(
        self,
        sample_size=100,
        initial_state=None,
        verbose=False,
    ):
        """
        Search for a MAP (maximum aposteriori) estimate by sampling (space exploration)
            This method returns a single-point estimate of the MAP of the distribution
        :param int sample_size:
        :param initial_state:
        :param bool verbose:
        """
        hmc_results = self.start_MCMC_sampling(
            sample_size=sample_size,
            initial_state=initial_state,
            full_diagnostics=full_diagnostics,
            verbose=verbose,
        )
        return hmc_results["map_estimate"]

    def generate_white_noise(self, size, truncate=True):
        """
        Generate a standard normal random vector of size `size` with values truncated
            at -/+3 if `truncate` is set to `True`

        :returns: a numpy array of size `size` sampled from a standard multivariate normal
            distribution of dimension `size` with mean 0 and covariance matrix equals
            an identity matrix.

        :remarks:
            - this function returns a numpy array of size `size` even if `size` is set to 1
        """
        # Sample (given the current internal random state), then reset the state, and truncate
        np_state = np.random.get_state()
        np.random.set_state(self.random_state)
        randn_vec = np.random.randn(size)
        self.random_state = np.random.get_state()
        np.random.set_state(np_state)

        if truncate:
            randn_vec[randn_vec > 3] = 3
            randn_vec[randn_vec < -3] = -3
        return randn_vec

    def mass_matrix_matvec(self, momentum):
        """
        Multiply the mass matrix (in the configurations) by the passed momentum
        """
        momentum = np.asarray(momentum).flatten()
        if momentum.size != self._CONFIGURATIONS["size"]:
            print(
                f"The passed momentum has invalid size;"
                f"received {momentum}, expected {self._CONFIGURATIONS['size']}"
            )
            raise TypeError

        return self._MASS_MATRIX.dot(momentum)

    def mass_matrix_inv_matvec(self, momentum):
        """
        Multiply the inverse of the mass matrix (in the configurations) by the passed momentum
        """
        momentum = np.asarray(momentum).flatten()
        if momentum.size != self._CONFIGURATIONS["size"]:
            print(
                f"The passed momentum has invalid size;"
                f"received {momentum}, expected {self._CONFIGURATIONS['size']}"
            )
            raise TypeError

        return self._MASS_MATRIX_INV.dot(momentum)

    def mass_matrix_sqrt_matvec(self, momentum):
        """
        Multiply the Square root (Lower Cholesky factor) of the mass matrix (in the configurations) by the passed momentum
        """
        momentum = np.asarray(momentum).flatten()
        if momentum.size != self._CONFIGURATIONS["size"]:
            print(
                f"The passed momentum has invalid size;"
                f"received {momentum}, expected {self._CONFIGURATIONS['size']}"
            )
            raise TypeError

        return self._MASS_MATRIX_SQRT.dot(momentum)

    def log_density(self, state):
        """
        Evaluate the value of the logarithm of the target unscaled posterior density function
        """
        val = self._LOG_DENSITY(state)
        try:
            val[0]
            val = np.asarray(val).flatten()[0]
        except:
            pass
        # if isinstance(val, np.ndarray) and val.size == 1: val = val.flatten()[0]
        return val

    def log_density_grad(self, state):
        """
        Evaluate the gradient of the logarithm of the target unscaled posterior density function
        """
        return self._LOG_DENSITY_GRAD(state)

    def potential_energy(self, state, verbose=False):
        """
        Evaluate the value of the potential energy at the given `state`
            The potential energy is the negative value of the logarithm of
            the unscaled posterior density function
        """
        if np.any(np.isnan(state)):
            if verbose:
                print("NaN values in the passed state")
                print(f"Received State:\n {repr(state)}")
            # raise ValueError
            return np.nan
        return -self.log_density(state)

    def potential_energy_grad(self, state, verbose=False):
        """
        Evaluate the gradient of the potential energy at the given `state`
            The potential energy is the negative value of the logarithm of
            the unscaled posterior density function
        """
        if np.any(np.isnan(state)):
            if verbose:
                print("NaN values in the passed state")
                print(f"Received State:\n {repr(state)}")
            # raise ValueError
            return np.nan
        return -self.log_density_grad(state)

    def kinetic_energy(self, momentum):
        """
        Evaluate the Kinetic energy of the posterior; this is independent from the state
            and is evaluated as the weighted l2 norm of the momentum
            (scaled by the inverse of hte mass matrix);
            This is half of the squared Mahalanobis distance of the Gaussian momentum

        :raises:
            - :py:class:`TypeError` is raised if the passed momentum has invalid shape/type/size
        """
        momentum = np.asarray(momentum).flatten()
        if momentum.size != self._CONFIGURATIONS["size"]:
            print(
                f"The passed momentum has invalid size;"
                f"received {momentum}, expected {self._CONFIGURATIONS['size']}"
            )
            raise TypeError

        return 0.5 * np.dot(momentum, self.mass_matrix_inv_matvec(momentum))

    def total_Hamiltonian(self, momentum, state):
        """
        Evaluate the value of the total energy function:
            Hamiltonian = kinetic energy + potential energy
        """
        return self.kinetic_energy(momentum) + self.potential_energy(state)

    def build_Hamiltonian_trajectory(
        self,
        momentum,
        state,
        step_size=_DEF_CONFIGURATIONS["symplectic_integrator_stepsize"],
        num_steps=_DEF_CONFIGURATIONS["symplectic_integrator_num_steps"],
        randomize_step_size=False,
    ):
        """
        Given the current momentum and state pair of the Hamiltonian system, generate a trajectory
            of (momentum, state).
        :param momentum:
        :param state:
        :param bool randomize_step_size: if `True` a tiny random number is added to the passed step
            size to help improve space exploration
        """
        # local copies
        momentum = np.asarray(momentum).flatten()
        state = np.asarray(state).flatten()

        trajectory = [(momentum, state)]
        # Loop over number of steps, for each step update current momentum and state then append to trajectory
        for _ in range(num_steps):
            trajectory.append(
                apply_symplectic_integration(
                    momentum=trajectory[-1][0],
                    state=trajectory[-1][1],
                )
            )
        return trajectory

    def apply_symplectic_integration(
        self,
        momentum,
        state,
        step_size=_DEF_CONFIGURATIONS["symplectic_integrator_stepsize"],
        num_steps=_DEF_CONFIGURATIONS["symplectic_integrator_num_steps"],
        randomize_step_size=False,
        symplectic_integrator="3-stage",
        SD=np.ones(50),
    ):
        """
        Apply one full step of size `step_size` of the symplectic integrator to the Hamiltonian system

        :parm momentum:
        :param state:
        :param int num_steps:
        :param float step_size:
        :param str symplectic_integrator: name of the symplectic integrator to use;
            acceptable are: 'verlet', 'leapfrog', '2-stage', '3-stage',
                where both 'leapfrog' and 'verlet' are equivalent
        :param bool randomize_step_size: if `True` a tiny random number is added to the passed step
            size to help improve space exploration

        :returns: a tuple (p, s) where p and s are the integrated (forward in time) momentum and state respectively
        """
        # local copies
        if np.any(np.isnan(momentum)):
            print(
                "Cannot apply symplectic integorator; NaN values found in the passed momentum"
            )
            raise ValueError
        if np.any(np.isnan(state)):
            print(
                "Cannot apply symplectic integorator; NaN values found in the passed state"
            )
            raise ValueError
        current_momentum = np.asarray(momentum).flatten()
        current_state = np.asarray(state).flatten()

        state_space_dimension = self._CONFIGURATIONS["size"]
        if not (current_momentum.size == current_state.size == state_space_dimension):
            print(
                f"The momentum and state must be of the same size as the underlying space dimnsion; "
                f"State size: {current_state.size}, "
                f"Momentum size: {current_momentum.size}, "
                f"Underlying space dimension: {state_space_dimension}"
            )
            raise TypeError

        # validate step size (and randomize if asked)
        if step_size <= 0:
            print(f"Step size of the symplectic integrator must be positive!!")
            raise ValueError

        if randomize_step_size:
            # random step size perturbation (update random state)
            np_state = np.random.get_state()
            np.random.set_state(self.random_state)
            u = (np.random.rand() - 0.5) * 0.4  # perturb step-size:
            self.random_state = np.random.get_state()
            np.random.set_state(np_state)

            h = (1 + u) * step_size
        else:
            h = step_size

        reflection_flag = False
        i = 0
        for _ in range(num_steps):
            #
            if re.match(r"\A(verlet|leapfrog)\Z", symplectic_integrator, re.IGNORECASE):
                # print(i, ": current state", current_state)
                # Update state
                proposed_state = current_state + (
                    0.5 * h
                ) * self.mass_matrix_inv_matvec(current_momentum)
                # print(i, ": proposed state", proposed_state)

                # Update momentum
                grad = self.potential_energy_grad(proposed_state)
                proposed_momentum = current_momentum - h * grad
                # print(i, ": proposed momentum", proposed_momentum)

                # Update state again
                proposed_state += (0.5 * h) * self.mass_matrix_inv_matvec(
                    proposed_momentum
                )
                # print(i, ": proposed state", proposed_state)
                i += 1

            elif re.match(
                r"\A2(-|_| )*stage(s)*\Z", symplectic_integrator, re.IGNORECASE
            ):
                a1 = 0.21132
                a2 = 1.0 - 2.0 * a1
                b1 = 0.5

                proposed_state = current_state + (a1 * h) * self.mass_matrix_inv_matvec(
                    current_momentum
                )

                grad = self.potential_energy_grad(proposed_state)
                proposed_momentum = current_momentum - (b1 * h) * grad

                proposed_state = proposed_state + (
                    a2 * h
                ) * self.mass_matrix_inv_matvec(proposed_momentum)

                grad = self.potential_energy_grad(proposed_state)
                proposed_momentum = proposed_momentum - (b1 * h) * grad

                proposed_state += (a1 * h) * self.mass_matrix_inv_matvec(
                    proposed_momentum
                )

            elif re.match(
                r"\A3(-|_| )*stage(s)*\Z", symplectic_integrator, re.IGNORECASE
            ):
                a1 = 0.11888010966548
                a2 = 0.5 - a1
                b1 = 0.29619504261126
                b2 = 1.0 - 2.0 * b1

                proposed_state = current_state + (a1 * h) * self.mass_matrix_inv_matvec(
                    current_momentum
                )

                grad = self.potential_energy_grad(proposed_state)
                proposed_momentum = current_momentum - (b1 * h) * grad

                proposed_state = proposed_state + (
                    a2 * h
                ) * self.mass_matrix_inv_matvec(proposed_momentum)

                grad = self.potential_energy_grad(proposed_state)
                proposed_momentum = proposed_momentum - (b2 * h) * grad

                proposed_state = proposed_state + (
                    a2 * h
                ) * self.mass_matrix_inv_matvec(proposed_momentum)

                grad = self.potential_energy_grad(proposed_state)
                proposed_momentum = proposed_momentum - (b1 * h) * grad

                proposed_state += (a1 * h) * self.mass_matrix_inv_matvec(
                    proposed_momentum
                )

            else:
                raise ValueError(
                    "Unsupported symplectic integrator %s" % symplectic_integrator
                )

            # Update current state and momentum

            if np.all(proposed_state > 0):
                current_momentum = proposed_momentum
                current_state = proposed_state
            else:
                negative_indices = np.where(proposed_state < 0)[0]
                # print("Indices of negative elements:", negative_indices)
                # print("before_change_proposed_state",proposed_state)
                proposed_state[negative_indices] = 0.01 * SD[negative_indices]
                current_momentum = proposed_momentum
                current_state = proposed_state
                reflection_flag = True
                print("reflection_flag", reflection_flag)
                # print("after_change_proposed_state",proposed_state)

        epsilon_indices = proposed_state <= 0.01 * SD

        return (proposed_momentum, proposed_state, epsilon_indices)

    def start_MCMC_sampling(
        self,
        sample_size,
        initial_state=None,
        randomize_step_size=False,
        full_diagnostics=False,
        verbose=False,
    ):
        """
        Start the HMC sampling procedure with initial state as passed.
        Use the underlying configurations for configuring the Hamiltonian trajectory, burn-in and mixin settings.

        :param int sample_size: number of smaple points to generate/collect from the predefined target distribution
        :param initial_state: initial point of the chain (any point that falls in the target distribution or near by it
            will result in faster convergence). You can try prior mean if this is used in a Bayesian approach
        :param bool randomize_step_size: if `True` a tiny random number is added to the passed step
            size to help improve space exploration
        :param bool full_diagnostics: if `True` all generated states will be tracked and kept for full disgnostics, otherwise,
            only collected samples are kept in memory
        :param bool verbose: screen verbosity
        """

        # Extract configurations from the configurations dictionary
        state_space_dimension = self._CONFIGURATIONS["size"]
        burn_in_steps = self._CONFIGURATIONS["burn_in"]
        mixing_steps = self._CONFIGURATIONS["mix_in"]
        symplectic_integrator = self._CONFIGURATIONS["symplectic_integrator"]
        hamiltonian_step_size = self._CONFIGURATIONS["symplectic_integrator_stepsize"]
        hamiltonian_num_steps = self._CONFIGURATIONS["symplectic_integrator_num_steps"]
        constraint_test = self._CONFIGURATIONS["constraint_test"]
        SD_theta = self._CONFIGURATIONS["SD_theta"]
        SD_gamma = self._CONFIGURATIONS["SD_gamma"]
        SD = np.concatenate((SD_theta, SD_gamma))
        liner, sliner = "=" * 53, "-" * 40
        if verbose:
            print("\n%s\nStarted Sampling\n%s\n" % (liner, liner))

        # Chain initial state
        if initial_state is None:
            initial_state = self.generate_white_noise(state_space_dimension)
        else:
            initial_state = np.array(initial_state).flatten()
            if initial_state.size != state_space_dimension:
                print(
                    f"Passed initial stae has invalid shape/size"
                    f"Passed initial state has size {initial_state.size}"
                    f"Expected size: {state_space_dimension}"
                )
                raise TypeError

        # Setup and construct the chain using HMC proposal:
        chain_length = burn_in_steps + sample_size * mixing_steps

        # print("initial_state",initial_state)
        # Initialize the chain
        current_state = initial_state.copy()  # initial state = ensemble mean

        # All generated sample points will be kept for testing and efficiency analysis
        chain_state_repository = [initial_state]
        proposals_repository = []
        acceptance_flags = []
        acceptance_probabilities = []
        uniform_random_numbers = []
        collected_ensemble = []
        map_estimate = None
        map_estimate_log_density = -np.infty

        count_positive = 0
        site_count_positive = np.zeros(initial_state.size)

        # Build the Markov chain
        start_time = time.time()  # start timing
        for chain_ind in range(chain_length):
            ## Proposal step :propose (momentum, state) pair
            # Generate a momentum proposal
            current_momentum = self.generate_white_noise(size=state_space_dimension)
            current_momentum = self.mass_matrix_sqrt_matvec(current_momentum)

            # print("current_momentum",current_momentum)
            # Advance the current state and momentum to propose a new pair:
            (
                proposed_momentum,
                proposed_state,
                epsilon_indices,
            ) = self.apply_symplectic_integration(
                momentum=current_momentum,
                state=current_state,
                num_steps=hamiltonian_num_steps,
                step_size=hamiltonian_step_size,
                randomize_step_size=randomize_step_size,
                symplectic_integrator=symplectic_integrator,
                SD=SD,
            )
            site_count_positive += epsilon_indices
            # if reflection_flag==True:
            #     count_positive+=1
            # print("site_count_positive",site_count_positive)

            #             print("proposed_momentum, proposed_state", proposed_momentum, proposed_state)

            ## MH step (Accept/Reject) proposed (momentum, state)
            # Calculate acceptance proabability
            # Total energy (Hamiltonian) of the extended pair (proposed_momentum,
            # Here, we evaluate the kernel of the posterior at both the current and the proposed state proposed_state)
            current_energy = self.total_Hamiltonian(
                momentum=current_momentum, state=current_state
            )
            constraint_violated = False

            # print("proposed_momentum",proposed_momentum,"current_momentum",current_momentum)
            # print("before_test_proposed_state",proposed_state,"before_test_current_state",current_state)

            if constraint_test is not None:
                if not constraint_test(proposed_state):
                    constraint_violated = True

            if constraint_violated:
                acceptance_probability = 0
                print("constraint violated")
            else:
                proposal_kinetic_energy = self.kinetic_energy(proposed_momentum)
                proposal_potential_energy = self.potential_energy(proposed_state)
                proposal_energy = proposal_kinetic_energy + proposal_potential_energy

                energy_loss = proposal_energy - current_energy

                # print("proposed_momentum",proposed_momentum,"current_momentum",current_momentum)
                # print("proposed_state",proposed_state,"current_state",current_state)
                # print("current_momentum",np.max(proposed_momentum),"current_momentum",np.max(current_momentum))
                # print("energy_loss",energy_loss,"proposal_energy",proposal_energy,"current_energy",current_energy,
                #         "proposal_kinetic_energy",proposal_kinetic_energy,"proposal_potential_energy",proposal_potential_energy)
                # print("current_kinetic_energy",self.kinetic_energy(current_momentum),"current_potential_energy",self.potential_energy(current_state))

                _loss_thresh = 1000
                if (
                    abs(energy_loss) >= _loss_thresh
                ):  # this should avoid overflow errors
                    if energy_loss < 0:
                        sign = -1
                    else:
                        sign = 1
                    energy_loss = sign * _loss_thresh
                    print("energy loss")
                acceptance_probability = np.exp(-energy_loss)

                acceptance_probability = min(acceptance_probability, 1.0)

                # Update Mode (Map Point Estimate)
                if -proposal_potential_energy > map_estimate_log_density:
                    map_estimate = proposed_state.copy()
                    map_estimate_log_density = -proposal_potential_energy

            # a uniform random number between 0 and 1
            np_state = np.random.get_state()
            np.random.set_state(self.random_state)
            uniform_probability = np.random.rand()
            self.random_state = np.random.get_state()
            np.random.set_state(np_state)

            # MH-rule
            if acceptance_probability > uniform_probability:
                current_state = proposed_state
                accept_proposal = True
            else:
                accept_proposal = False

            if verbose:
                print(
                    f"\rHMC Iteration [{chain_ind+1:4d}/{chain_length:4d}]; Accept Prob: {acceptance_probability:3.2f}; --> Accepted? {accept_proposal}",
                    end="  ",
                )

            #
            if chain_ind >= burn_in_steps and chain_ind % mixing_steps == 0:
                collected_ensemble.append(current_state.copy())

            # Update Results Repositories:
            if full_diagnostics:
                proposals_repository.append(proposed_state)
                acceptance_probabilities.append(acceptance_probability)
                uniform_random_numbers.append(uniform_probability)
                #
                if accept_proposal:
                    acceptance_flags.append(1)
                else:
                    acceptance_flags.append(0)
                chain_state_repository.append(np.squeeze(current_state))

        print("site_count_positive", site_count_positive)
        # Stop timing
        chain_time = time.time() - start_time

        # ------------------------------------------------------------------------------------------------

        # Now output diagnostics and show some plots :)
        if full_diagnostics:
            chain_diagnostics = self.mcmc_chain_diagnostic_statistics(
                proposals_repository=proposals_repository,
                chain_state_repository=chain_state_repository,
                collected_ensemble=collected_ensemble,
                acceptance_probabilities=acceptance_probabilities,
                uniform_probabilities=uniform_random_numbers,
                acceptance_flags=acceptance_flags,
                map_estimate=map_estimate,
            )
        else:
            chain_diagnostics = None

        #
        # ======================================================================================================== #
        #                Output sampling diagnostics and plot the results for 1 and 2 dimensions                   #
        # ======================================================================================================== #
        #
        if verbose:
            print("MCMC sampler:")
            print(f"Time Elapsed for MCMC sampling: {chain_time} seconds")
            if chain_diagnostics is not None:
                print(f"Acceptance Rate: {chain_diagnostics['acceptance_rate']:.2f}")

        sampling_results = dict(
            chain_state_repository=chain_state_repository,
            collected_ensemble=collected_ensemble,
            proposals_repository=proposals_repository,
            acceptance_flags=acceptance_flags,
            acceptance_probabilities=acceptance_probabilities,
            uniform_random_numbers=uniform_random_numbers,
            chain_diagnostics=chain_diagnostics,
            map_estimate=map_estimate,
            map_estimate_log_density=map_estimate_log_density,
            chain_time=chain_time,
            count_positive=count_positive,
            site_count_positive=site_count_positive,
        )
        return sampling_results

    def _mcmc_chain_diagnostic_statistics(
        self,
        proposals_repository,
        chain_state_repository,
        uniform_probabilities,
        acceptance_probabilities,
        collected_ensemble,
        acceptance_flags=None,
    ):
        """
        Return diagnostic statistics of the chain such as the rejection rate, acceptance ratio, etc.
        """
        if acceptance_flags is None:
            acceptance_flags = np.asarray(
                acceptance_probabilities >= uniform_probabilities, dtype=np.int
            )
        else:
            acceptance_flags = np.asarray(acceptance_flags)
        acceptance_rate = (
            float(acceptance_flags.sum()) / np.size(acceptance_flags) * 100.0
        )
        rejection_rate = 100.0 - acceptance_rate

        # TODO: Add More:
        # effective sample size, autocorrelation, etc.

        # Return all diagonistics in a dictionary
        chain_diagnositics = dict(
            acceptance_rate=acceptance_rate,
            rejection_rate=rejection_rate,
        )
        return chain_diagnositics

    @property
    def random_state(self):
        """Get a handle of the current internal random state"""
        return self._RANDOM_STATE

    @random_state.setter
    def random_state(self, value):
        """Update the internal random state"""
        try:
            np_state = np.random.get_state()
            np.random.set_state(value)
            self._RANDOM_STATE = value
            np.random.set_state(np_state)
        except:
            print("Invalid random state passed of type '{0}'".format(type(value)))
            raise TypeError


## Simple interfaces (to generate instances from classes developed here).
def create_hmc_sampler(
    size,
    log_density,
    log_density_grad=None,
    burn_in=100,
    mix_in=10,
    symplectic_integrator="verlet",
    symplectic_integrator_stepsize=1e-2,
    symplectic_integrator_num_steps=20,
    #    mass_matrix=(5000,1),
    random_seed=1011,
    constraint_test=None,
    mass_matrix_gamma=np.ones(25),
    mass_matrix_theta=np.ones(25),
    mass_matrix_gamma_scale=1.0,
    mass_matrix_theta_scale=1.0,
):
    """
    Given the size of the target space, and a function to evalute log density,
        create and return an :py:class:`HMCSampler` instance/object to generate samples using HMC sampling approach.
        Configurations/settings can be updated after inistantiation

    This function shows how to create :py:class:`HMCSampler` instances (with some or all configurations passed)
    """

    configs = dict(
        size=size,
        log_density=log_density,
        log_density_grad=log_density_grad,
        burn_in=burn_in,
        mix_in=mix_in,
        symplectic_integrator=symplectic_integrator,
        symplectic_integrator_stepsize=symplectic_integrator_stepsize,
        symplectic_integrator_num_steps=symplectic_integrator_num_steps,
        # mass_matrix=mass_matrix,
        random_seed=random_seed,
        constraint_test=constraint_test,
        SD_gamma=np.sqrt(1 / mass_matrix_gamma),
        SD_theta=np.sqrt(1 / mass_matrix_theta),
        mass_matrix_gamma=mass_matrix_gamma * mass_matrix_gamma_scale,
        mass_matrix_theta=mass_matrix_theta * mass_matrix_theta_scale,
    )
    return HMCSampler(configs)


if __name__ == "__main__":
    # This must exist for multiprocessing to work
    pass
