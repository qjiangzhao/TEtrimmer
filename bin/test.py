import os
import unittest
from TE_Trimmer import main  # Import your main function
from click.testing import CliRunner
import click


class TestMainFunction(unittest.TestCase):

    def test_main_function(self):
        # Set up test data or input arguments
        input_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tests/test_input.fa')
        genome_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tests/test_genome.fasta')
        output_dir = 'test_output_dir'
        os.makedirs(output_dir, exist_ok=True)
        if os.listdir(output_dir):
            click.echo(
                f"WARNING: The output directory {output_dir} is not empty. Please empty your output directory or "
                f"choose another empty directory\n")
            # Stop the whole program when the output directory is not empty
            return
        
        args = [
            '--input_file', input_file,
            '--genome_file', genome_file,
            '--output_dir', output_dir,
            '--species', 'fungi',
            '--classify_unknown'
        ]
        runner = CliRunner()
        result = runner.invoke(main, args)

        # Add assertions to check the results
        # check input fasta
        files_in_folder = os.listdir(f"{output_dir}/Single_fasta_files")
        expected_file_count = 8
        self.assertEqual(len(files_in_folder), expected_file_count, f"Expected {expected_file_count} files in the {output_dir}/Single_fasta_files, but found {len(files_in_folder)} files.")

        # check summary.txt
        self.assertTrue(os.path.exists('test_output_dir/summary.txt'), 'summary.txt not found')
        summary_size = os.path.getsize(os.path.join(output_dir, 'summary.txt'))
        self.assertGreater(summary_size, 800, 'summary.txt is not complete')

        # check TE_Trimmer_consensus.fasta
        self.assertTrue(os.path.exists('test_output_dir/TE_Trimmer_consensus.fasta'), 'TE_Trimmer_consensus.fasta not found')
        te_trimmer_size = os.path.getsize(os.path.join(output_dir, 'TE_Trimmer_consensus.fasta'))
        self.assertGreater(te_trimmer_size, 10000, 'TE_Trimmer_consensus.fasta is not complete')

        # check proof_annotation
        contents = os.listdir(f"{output_dir}/TE_Trimmer_for_proof_annotation/Need_check_annotation")
        self.assertGreater(len(contents), 0, f"{output_dir}/TE_Trimmer_for_proof_annotation/Need_check_annotation is empty.")


if __name__ == '__main__':
    unittest.main()
