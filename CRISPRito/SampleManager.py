import pandas as pd
import uuid
import os

class SampleManager:
    def __init__(self, input_file, output_dir="sample_outputs"):
        """
        Initializes the SampleManager with an input file and output directory.
        
        :param input_file: Path to the input file containing sample data.
        :param output_dir: Directory where individual sample files will be saved.
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.samples = None  # Will store DataFrame after loading
    
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

    def load_samples(self):
        """Reads the input file into a pandas DataFrame and assigns unique IDs."""
        self.samples = pd.read_csv(self.input_file, sep="\t")
        
        # Assign unique identifiers
        self.samples["id"] = [str(uuid.uuid4()) for _ in range(len(self.samples))]

    def write_samples(self):
        """Writes each row to a separate text file using its unique ID as the filename."""
        if self.samples is None:
            raise ValueError("Samples must be loaded first using load_samples().")

        for _, row in self.samples.iterrows():
            file_path = os.path.join(self.output_dir, f"{row['id']}.txt")
            
            with open(file_path, "w") as f:
                f.write(row.to_json(indent=4))  # Write row data as JSON for readability

    def process(self):
        """Runs the full pipeline: loading samples and writing output files."""
        self.load_samples()
        self.write_samples()

