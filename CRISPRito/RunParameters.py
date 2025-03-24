import os
import shutil

class RunParameters:
    def __init__(self, output_dir, overwrite=False):
        """
        Initializes RunParameters and manages the output directory.

        :param output_dir: Path to the output directory.
        :param overwrite: If True, deletes and recreates the directory if it exists.
        """
        self.output_dir = output_dir
        self.overwrite = overwrite
        self._handle_output_dir()

    def _handle_output_dir(self):
        """Creates or resets the output directory based on overwrite setting."""
        if os.path.exists(self.output_dir):
            if self.overwrite:
                shutil.rmtree(self.output_dir)  # Delete existing directory
                print(f"Output directory '{self.output_dir}' exists. Overwriting...")
                os.makedirs(self.output_dir)
            else:
                print(f"Output directory '{self.output_dir}' already exists. Keeping existing files.")
        else:
            os.makedirs(self.output_dir)
            print(f"Created output directory: {self.output_dir}")

    def get_output_dir(self):
        """Returns the output directory path."""
        return self.output_dir
