import unittest
import tempfile
import py_compile
import subprocess
import shutil
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
        cls.input_file = os.path.join(cls.test_data_dir, 'Test.vaccine.results.input.fa')
        cls.method = 'ann'
        cls.keep_tmp = 'True'
        cls.allele = 'H-2-Kb'
        cls.epitope_length = str(8)
        cls.seed = 'True'

    # def test_vaccine_design_compiles(self):
    #     self.assertTrue(py_compile.compile(self.executable))

    def test_vaccine_design_runs_and_produces_expected_output(self):
        output_dir = os.path.join(self.test_data_dir, tempfile.mkdtemp())
        tempfile.tempdir = os.path.join(self.test_data_dir, output_dir)

        subprocess.call([self.python,
                           self.executable,
                           self.test_run_name,
                           self.input_file,
                           self.method,
                           self.allele,
                           '-o', tempfile.gettempdir(),
                           '-e', self.epitope_length,
                           '-k', self.keep_tmp,
                           '-s', self.seed], shell=False)

        self.assertTrue(cmp(
            os.path.join(tempfile.gettempdir(), self.test_run_name, self.test_run_name + '_results.fa'),
            os.path.join(self.test_data_dir, "Test.vaccine.results.output.fa")
        ))

        self.assertTrue(cmp(
            os.path.join(tempfile.gettempdir(), self.test_run_name, 'tmp', self.test_run_name + '_epitopes.fa'),
            os.path.join(self.test_data_temp_dir, 'Test.vaccine.design.epitopes.comb.fa')
        ))

        self.assertTrue(cmp(
            os.path.join(tempfile.gettempdir(), self.test_run_name, 'tmp', self.test_run_name + '_iedb_out.csv'),
            os.path.join(self.test_data_temp_dir, 'Test.vaccine.design.iedb.results.csv')
        ))

        shutil.rmtree(output_dir)
        
