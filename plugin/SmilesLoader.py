import nanome
from nanome.api.structure import Complex
from nanome.api.ui import Menu
from nanome.util.enums import Integrations, NotificationTypes

from rdkit import Chem
from rdkit.Chem import AllChem

import tempfile


class SmilesLoader(nanome.PluginInstance):
    def start(self):
        self.integration.export_smiles = self.integration_export
        self.integration.import_smiles = self.integration_import

        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)

        menu = Menu()
        self.menu = menu

        menu.title = 'SMILES'
        menu.width = 0.5
        menu.height = 0.2

        menu.root.padding_type = menu.root.PaddingTypes.ratio
        menu.root.set_padding(top=0.05, down=0.02, left=0.01, right=0.01)

        ln_inp = menu.root.create_child_node()
        ln_inp.forward_dist = 0.001
        inp_smiles = ln_inp.add_new_text_input()
        inp_smiles.text_size = 0.4
        inp_smiles.max_length = 1000
        inp_smiles.placeholder_text = 'Enter SMILES'
        self.inp_smiles = inp_smiles

        ln_btn = menu.root.create_child_node()
        ln_btn.forward_dist = 0.001
        btn_load = ln_btn.add_new_button('Load')
        btn_load.text.value.unusable = 'Loading...'
        btn_load.disable_on_press = True
        btn_load.register_pressed_callback(self.load_from_input)

    def on_run(self):
        self.menu.enabled = True
        self.update_menu(self.menu)

    def integration_export(self, request):
        complexes = request.get_args()
        strings = self.complexes_to_strings(complexes)
        request.send_response(strings)

    def integration_import(self, request):
        strings = request.get_args()
        complexes = []

        for string in strings:
            mols = []
            with Chem.SDWriter(self.temp_sdf.name) as w:
                for smiles in string.split('\n'):
                    mol = self.smiles_to_mol(smiles)
                    if mol is None:
                        continue
                    mols.append(mol)
                    w.write(mol)

            if not mols:
                continue

            complex = Complex.io.from_sdf(path=self.temp_sdf.name)
            complex.name = 'SMILES ' + AllChem.CalcMolFormula(mols[0])
            complexes.append(complex)

        request.send_response(complexes)

    def load_from_input(self, btn=None):
        smiles = self.inp_smiles.input_text
        mol = self.smiles_to_mol(smiles)

        if mol is None:
            msg = 'RDKit could not parse the SMILES.'
            self.send_notification(NotificationTypes.error, msg)
        else:
            with Chem.SDWriter(self.temp_sdf.name) as w:
                w.write(mol)
            name = 'SMILES ' + AllChem.CalcMolFormula(mol)
            self.send_files_to_load((self.temp_sdf.name, name))

        self.inp_smiles.input_text = ''
        btn.unusable = False
        self.update_content([self.inp_smiles, btn])

    def smiles_to_mol(self, smiles):
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return None

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0)
        return mol

    def complexes_to_strings(self, complexes):
        strings = []

        for complex in complexes:
            complex.io.to_sdf(self.temp_sdf.name)
            lines = []

            for mol in Chem.SDMolSupplier(self.temp_sdf.name):
                if mol is None:
                    continue
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                lines.append(smiles)

            strings.append('\n'.join(lines))

        return strings


def main():
    plugin = nanome.Plugin('SMILES Loader', 'A Nanome Plugin to load from SMILES string using RDKit', 'Files', False, integrations=[Integrations.smiles])
    plugin.set_plugin_class(SmilesLoader)
    plugin.run()


if __name__ == '__main__':
    main()
