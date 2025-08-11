import re


def safe_col(col):
    return re.sub("[ ,\:\"'\`\^\#\[\]\{\}\<\>\?\!\*\&\$\%\\\/\|]", "_", col)
