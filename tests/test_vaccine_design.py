import unittest
import tempfile
import py_compile
import subprocess
from filecmp import cmp
import os
import sys

#python -m unittest tests/test_vaccine_design.py
class TestVaccineDesign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.python = sys.executable
        cls.executable_dir = os.path.join(base_dir, 'pvacseq', 'lib')
        cls.executable = os.path.join(cls.executable_dir, 'vaccine_design.py')
        cls.test_run_name = 'test_vaccine_design_produces_expected_output'
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data', 'vaccine_design')
        cls.test_data_temp_dir = os.path.join(base_dir, 'tests', 'test_data', 'vaccine_design', 'tmp')
        cls.test_tmp_epitope_file = os.path.join(cls.test_data_temp_dir, cls.test_run_name + '_epitopes.fa')
        cls.test_tmp_iedb_out_file = os.path.join(cls.test_data_temp_dir, cls.test_run_name + '_iedb_out.csv' )
        cls.input_file = os.path.join(cls.test_data_dir, 'peptides.fa')
        cls.method = 'ann'
        cls.keep_tmp = 'True'
        cls.allele = 'H-2-Kb'
        cls.epitope_length = str(8)

    def test_vaccine_design_compiles(self):
        compiled_script_path = py_compile.compile(self.executable)
        self.assertTrue(compiled_script_path)

    def test_vaccine_design_runs_and_produces_expected_output(self):
        output_dir = tempfile.mkdtemp()
        if self.python is None:
            raise Exception("Python could not be found!")
        try:
            subprocess.check_call([self.python,
                  self.executable,
                  self.test_run_name,
                  self.input_file,
                  self.method,
                  self.allele,
                  '-o ' + output_dir,
                  '-e ' + self.epitope_length,
                  '-k ' + self.keep_tmp])
        except subprocess.CalledProcessError as exc:
            print ("error code", exc.returncode, exc.output)

        self.assertTrue(cmp(
            os.path.join(output_dir, self.test_run_name + '_results.fa'),
            os.path.join(self.test_data_dir, "Test.vaccine.results.output.fa")
        ))

        self.assertTrue(cmp(
            os.path.join(output_dir, 'tmp', self.test_run_name + '_epitopes.fa'),
            os.path.join(self.test_data_temp_dir, 'Test.vaccine.design.epitopes.comb.fa')
        ))

        self.assertTrue(cmp(
            os.path.join(output_dir, 'tmp', self.test_run_name + '_iedb_out.csv'),
            os.path.join(self.test_data_temp_dir, 'Test.vaccine.design.iedb.results.csv')
        ))

        # Not sure we need to compare log files. Since simulated annealing uses random values, its likely
        #   that some values will not exactly match between runs
        self.assertAlmostEqual(
            os.path.join( output_dir, self.test_run_name + '.log'),
            os.path.join(self.test_data_temp_dir, 'Test.vaccine.design.log')
        )

