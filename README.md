# Project Amazon

## Installation

0. Clone git repository and move into `project-amazon/`
1. Create a new virtual environment
```
python -m venv venv
```
2. Install dependencies
```
python -m pip install -e '.[notebooks, dev]'
```

## Contributing
0. Open a new git branch
```
git checkout -b <new branch name>
```
1. Create new code changes

2. Stage changed files
```
git add <names of changed files>
```

3. Commit changed files
```
git commit
```
4. After several commits, push commits to remote (if it is the first time pushing this branch add the --set-upstrea)
```
git push
```

5. Submit a pull request on GitHub