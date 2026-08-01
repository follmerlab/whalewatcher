import unittest

from orca_vib_viewer import _short_filename


class ShortFilenameTests(unittest.TestCase):
    def test_windows_path(self):
        self.assertEqual(
            _short_filename(r"C:\Users\chemist\modes\sample.out"),
            "sample.out",
        )

    def test_posix_path(self):
        self.assertEqual(
            _short_filename("/home/chemist/modes/sample.out"), "sample.out"
        )


if __name__ == "__main__":
    unittest.main()
