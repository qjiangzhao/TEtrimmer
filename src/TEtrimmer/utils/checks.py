import logging
import shutil
import sys


def check_tools(required_tools=[], optional_tools=[]):
    """
    Check if required and optional tools are available on the system's PATH.

    Args:
        required_tools (list): List of required tool names.
        optional_tools (list): List of optional tool names.

    Raises:
        RuntimeError: If any required tool is not found.
    """
    missing_required_tools = []

    def print_message(tool, path, color_code):
        """
        Print a message to stderr with the tool name and path in the specified color.

        Args:
            tool (str): The name of the tool.
            path (str): The path to the tool.
            color_code (str): The ANSI color code for the message.
        """
        tool_padded = tool.ljust(15)
        if path:
            message = f'{color_code}{tool_padded}\t{path}\033[0m'
        else:
            message = f'{color_code}{tool_padded}\tNOT FOUND\033[0m'
        print(message, file=sys.stderr)

    # Check required tools
    print('Checking for dependencies:')
    for tool in required_tools:
        path = shutil.which(tool)
        if path:
            print_message(tool, path, '\033[92m')  # Green
        else:
            print_message(tool, None, '\033[91m')  # Red
            missing_required_tools.append(tool)

    # Check optional tools
    for tool in optional_tools:
        path = shutil.which(tool)
        if path:
            print_message(tool, path, '\033[92m')  # Green
        else:
            print_message(tool, None, '\033[93m')  # Yellow

    # Raise error if any required tool is missing
    if missing_required_tools:
        error_message = 'ERROR: Some required tools could not be found: ' + ', '.join(
            missing_required_tools
        )
        logging.error(error_message)
        raise RuntimeError(
            'Missing required tools: ' + ', '.join(missing_required_tools)
        )
