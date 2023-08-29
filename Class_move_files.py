import os
import re
import shutil

class FileCopier:
    def __init__(self, source_dir, target_dir, patterns):
        """
        Based on regular expression pattern, move files from one folder to the other
        :param source_dir: str, directory contains all files to copy
        :param target_dir: str, absolute path to store copy files
        :param patterns: list, regular expression patterns
        """
        self.source_dir = source_dir
        self.target_dir = target_dir
        self.patterns = patterns

    def move_files(self):
        # Make target_dir if it doesn't exist
        if not os.path.exists(self.target_dir):
            os.makedirs(self.target_dir)

        files = os.listdir(self.source_dir)
        for file in files:
            for pattern in self.patterns:
                if re.match(pattern, file):
                    shutil.copy(os.path.join(self.source_dir, file), self.target_dir)

# Example usage:
# patterns = [r'\w+\.txt', r'\w+\.jpg']
# mover = FileMover("source_directory", "target_directory", patterns)
# mover.move_files()
