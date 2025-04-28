

def process_samples():



if __name__ == '__main__':
	import argparse


    parser = argparse.ArgumentParser(description="Process CRISPRito samples.")
    parser.add_argument("samplesheet", help="Path to the sample sheet CSV.")
    parser.add_argument("--output_dir", default="CRISPRito_output", help="Directory to save output.")

    args = parser.parse_args()
    process_samples(args.samplesheet, args.output_dir)

 	process_samples(args.samplesheet, args.output_dir)	