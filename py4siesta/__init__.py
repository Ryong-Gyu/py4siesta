"""py4siesta package."""

def main(*args, **kwargs):
    from .cli import main as cli_main

    return cli_main(*args, **kwargs)

__all__ = ["main"]
