import nanome
from nanome.api.ui import Menu
from nanome.util import enums

from rdkit import Chem
from rdkit.Chem import AllChem

import tempfile


class SmilesLoader(nanome.PluginInstance):
    def start(self):
        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf')

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
        btn = ln_btn.add_new_button('Load')
        btn.register_pressed_callback(self.load_smiles)

    def on_run(self):
        self.menu.enabled = True
        self.update_menu(self.menu)

    def load_smiles(self, btn=None):
        smiles = self.inp_smiles.input_text.strip()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            self.send_notification(enums.NotificationTypes.error, 'RDKit was unable to load this SMILES')
            return

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0)

        with open(self.temp_sdf.name, 'w') as f:
            f.write(Chem.MolToMolBlock(mol))

        self.send_files_to_load((self.temp_sdf.name, smiles.replace('/', '_')))
        self.inp_smiles.input_text = ''
        self.update_content(self.inp_smiles)


def main():
    plugin = nanome.Plugin('SMILES Loader', 'A Nanome Plugin to load from SMILES string using RDKit', 'Files', False)
    plugin.set_plugin_class(SmilesLoader)
    plugin.run()


if __name__ == '__main__':
    main()
