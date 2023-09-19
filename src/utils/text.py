class TextColor:
    """Class with text colors and font settings for printing"""

    PURPLE = "\033[95m"
    CYAN = "\033[96m"
    DARKCYAN = "\033[36m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


def decorate_text(text):
    sep = "*" * 48
    return f"""{TextColor.BOLD}{TextColor.DARKCYAN}\n\n
    {sep}\n
    \t{text}\n
    {sep}\n
    {TextColor.END}\n
    """
