#
# GAMS - General Algebraic Modeling System Python API
#
# Copyright (c) 2017-2024 GAMS Development Corp. <support@gams.com>
# Copyright (c) 2017-2024 GAMS Software GmbH <support@gams.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

from gams.tools.toolcollection.tooltemplate import ToolTemplate
from gams.connect import ConnectDatabase
import gams.transfer as gt


class Csvwrite(ToolTemplate):

    def __init__(self, system_directory, tool):
        super().__init__(system_directory, tool)
        self.title = "csvwrite: This tool exports a GAMS symbol to a CSV file."
        self.add_namedargdef(
            "allFields=<boolean>",
            "str",
            "Specify whether all the attributes (level, marginal, lower, upper, and scale) of a variable or an equation are written to the CSV. By default only the level will be written.",
            argdefault=False,
        )
        self.add_namedargdef(
            "decimalSep=<string>",
            "str",
            "Specify a decimal separator.",
            argdefault="period",
        )
        self.add_namedargdef(
            "dFormat=<string>",
            "str",
            "Specify the numerical format in the output file.",
            argdefault="normal",
        )
        self.add_namedargdef(
            "fieldSep=<string>", 
            "str", 
            "Specify a field separator.", 
            argdefault="comma",
        )
        self.add_namedargdef(
            "file=<string>",
            "fnWriteable",
            "Specify the name for the CSV file.",
        )
        self.add_namedargdef(
            "gdxIn=<gdx_filename>",
            "fnExist",
            "Specify the input GDX file",
            shell_req=True,
        )
        self.add_namedargdef(
            "header=<boolean>",
            "str",
            "Specify the header used as the column names.",
            argdefault=True,
        )
        self.add_namedargdef(
            "id=<string>",
            "str",
            "Specify a symbol name for the Connect database.",
        )
        self.add_namedargdef(
            "quoting=<int>",
            "int",
            "Control field quoting behavior.",
            argdefault=0
        )
        self.add_namedargdef(
            "setHeader=<string>",
            "str",
            "Specify a string that will be used as the header.",
            argdefault=None,
        )
        self.add_namedargdef(
            "skipText=<boolean>",
            "str",
            "Specify if the set element text is written.",
            argdefault=False,
        )
        self.add_namedargdef(
            "unstack=<boolean>",
            "str",
            "Specify if the last dimension will be unstacked to the header row.",
            argdefault=False,
        )

    @staticmethod
    def check_duplicate_in_named_args(tuple_list):
        seen_elements = set()
        duplicates = set()

        for tup in tuple_list:
            if tup[0] in seen_elements:
                duplicates.add(tup[0])
            else:
                seen_elements.add(tup[0])

        return ", ".join(duplicates)

    def check_bool_args(self, key):
        """
        Helper function to convert the argVal of a boolean type argument to Boolean True/False.
        Return default value if not specified.

        Raise Exception if the input is not y/n.
        """
        key = key.lower()
        if key in self.namedargs:
            value = self.namedargs_val(key)
            if value.lower() in ["y", "yes"]:
                return True
            elif value.lower() in ["n", "no"]:
                return False
            self.tool_error(f"Wrong flag: {key}={value}", print_help=False)

        return self.namedargs_val(key)  # return default

    def get_default_sep_values(self, key):
        separators = {
            "fieldsep": {"comma": ",", "semicolon": ";", "tab": "\t"},
            "decimalsep": {"period": ".", "comma": ","},
        }
        default_sep = {"fieldsep": ",", "decimalsep": "."}

        if key in self.namedargs:
            try:
                return separators[key][self.namedargs_val(key)]
            except:
                self.tool_error(
                    f"Wrong {key} input: {self.namedargs_val(key)}", print_help=False
                )

        return default_sep[key]

    def _execute_projection(
        self, cdb: ConnectDatabase, sym_name: str, all_fields: bool
    ):
        dim = cdb.container[sym_name].dimension
        if dim == 0:

            cdb.execute(
                {
                    "Projection": {
                        "name": f"{sym_name}.all" if all_fields else f"{sym_name}.l",
                        "newName": f"{sym_name}_all" if all_fields else f"{sym_name}_l",
                    }
                }
            )
        else:
            dom = ",".join(f"d{i}" for i in range(dim))
            cdb.execute(
                {
                    "Projection": {
                        "name": (
                            f"{sym_name}.all({dom})"
                            if all_fields
                            else f"{sym_name}.l({dom})"
                        ),
                        "newName": (
                            f"{sym_name}_all({dom})"
                            if all_fields
                            else f"{sym_name}_l({dom})"
                        ),
                    }
                }
            )
        if all_fields:
            return {"name": f"{sym_name}_all", "unstack": True}
        else:
            return {"name": f"{sym_name}_l"}

    def execute(self):
        if self.dohelp():
            return

        self.process_args()

        check_dupe_opts = self.check_duplicate_in_named_args(self.namedargs_list)
        if check_dupe_opts:
            self.tool_error(f"Duplicate option: {check_dupe_opts}", print_help=False)

        if not self.namedargs_val("id"):
            self.tool_error("Parameter >id< not set.")
        else:
            sym_name = self.namedargs_val("id")

        if self.namedargs_val("file"):
            csv_file = self.namedargs_val("file")
        elif not self.namedargs_val("file") and self.namedargs_val("gdxin"):
            csv_file = self.namedargs_val("gdxin").rsplit(".gdx", 1)[0] + ".csv"
        else:
            self.tool_error(f"Option >file< not specified.")

        cdb = ConnectDatabase(self._tools._system_directory, ecdb=self._tools._ecdb)
        self.read_id_inputs(cdb.container, inputs=sym_name)

        if self.namedargs_val("dformat") != "normal":
            _gdx_file = self.namedargs_val("gdxin")
            if not _gdx_file:
                import tempfile

                with tempfile.NamedTemporaryFile(
                    mode="w", dir=".", delete=False, suffix=".gdx"
                ) as fp:
                    cdb.container.write(fp.name, sym_name)
                    _gdx_file = fp.name

            from subprocess import run, PIPE
            from shutil import which

            if which("gdxdump"):
                map_csvwriter_to_gdxdump_args = {
                    "setheader": "header",
                    "unstack": "cdim",
                    "fieldsep": "delim",
                    "dformat": "dformat",
                    "decimalsep": "decimalsep",
                }
                cmd = [
                    "gdxdump",
                    _gdx_file,
                    f"output={csv_file}",
                    "format=csv",
                    f"symb={sym_name}",
                ]
                if self.check_bool_args("allfields"):
                    cmd.append("CSVAllFields")

                if not self.check_bool_args("skiptext"):
                    cmd.append("CSVSetText")

                if not self.check_bool_args("header"):
                    cmd.append("noHeader")

                for key in self.namedargs_list:
                    csvwriter_args = key[0]
                    if csvwriter_args in map_csvwriter_to_gdxdump_args:
                        cmd.append(
                            f"{map_csvwriter_to_gdxdump_args[csvwriter_args]}={self.namedargs_val(csvwriter_args)}"
                        )
                cmd_res = run(
                    cmd,
                    stdout=PIPE,
                    stderr=PIPE,
                    universal_newlines=True,
                )
                if not self.namedargs_val("gdxin"):
                    import os

                    os.remove(_gdx_file)
                if cmd_res.returncode != 0:
                    self.tool_error(
                        f"Error occured while running GDXDUMP utility. Error: {cmd_res.stdout}"
                    )
            else:
                self.tool_error(
                    f"GDXDUMP utility not found on the system.", print_help=False
                )
        else:
            csvwriter_params = {
                "file": csv_file,
                "name": sym_name,
                "trace": self.namedargs_val("trace"),
                "header": self.check_bool_args("header"),
                "unstack": self.check_bool_args("unstack"),
                "decimalSeparator": self.get_default_sep_values(key="decimalsep"),
                "fieldSeparator": self.get_default_sep_values(key="fieldsep"),
                "quoting": self.namedargs_val("quoting"),
                "setHeader": self.namedargs_val("setheader"),
                "skipText": self.check_bool_args("skiptext"),
            }
            if (
                csvwriter_params["decimalSeparator"]
                == csvwriter_params["fieldSeparator"]
            ):
                # GDXDUMP does not allow same separators
                self.tool_error(
                    f"fieldSep and decimalSep characters should be different."
                )
            if isinstance(cdb.container[sym_name], gt.Variable) or isinstance(
                cdb.container[sym_name], gt.Equation
            ):
                csvwriter_params.update(
                    self._execute_projection(
                        cdb, sym_name, self.check_bool_args("allFields")
                    )
                )
            try:
                cdb.execute({"CSVWriter": csvwriter_params})
            except Exception as e:
                self.tool_error(f"{e.__class__.__name__} : {e}", print_help=False)
