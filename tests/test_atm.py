import shutil
import os
import yaml

curr_dir = os.path.dirname(os.path.abspath(__file__))


def _test_abfe_structprep(tmp_path):
    from atom_openmm.abfe_structprep import abfe_structprep

    run_dir = os.path.join(tmp_path, "1STP")
    shutil.copytree(os.path.join(curr_dir, "1STP"), run_dir)

    abfe_structprep(os.path.join(run_dir, "1STP_asyncre.yaml"))


def _test_abfe_production(tmp_path):
    from atom_openmm.abfe_production import abfe_production

    run_dir = os.path.join(tmp_path, "1STP_equil")
    shutil.copytree(os.path.join(curr_dir, "1STP_equil"), run_dir)

    abfe_production(os.path.join(run_dir, "1STP_asyncre.yaml"))


def _test_rbfe_structprep(tmp_path):
    from atom_openmm.rbfe_structprep import rbfe_structprep
    from openmm import unit
    from openmm import app
    from openmm.app.amberprmtopfile import AmberPrmtopFile
    from openmm import XmlSerializer
    import numpy as np

    run_dir = os.path.join(tmp_path, "QB_A08_A07")
    shutil.copytree(os.path.join(curr_dir, "QB_A08_A07"), run_dir)

    box_vectors = np.array([[80.1356, 0, 0], [0, 86.0443, 0], [0, 0, 81.4926]])
    a = unit.Quantity(box_vectors[0] * unit.angstrom)
    b = unit.Quantity(box_vectors[1] * unit.angstrom)
    c = unit.Quantity(box_vectors[2] * unit.angstrom)
    box_vectors = (a, b, c)
    omm_structure = AmberPrmtopFile(
        os.path.join(run_dir, "QB_A08_A07.prmtop"), periodicBoxVectors=box_vectors
    )
    system = omm_structure.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=9 * unit.angstrom,
        constraints=app.HBonds,
        hydrogenMass=1.5 * unit.amu,
        rigidWater=True,
        switchDistance=7 * unit.angstrom,
    )
    with open(os.path.join(run_dir, "QB_A08_A07_sys.xml"), "w") as f:
        f.write(XmlSerializer.serialize(system))

    rbfe_structprep(os.path.join(run_dir, "QB_A08_A07_asyncre.yaml"))


def _test_rbfe_structprep_from_keywords(tmp_path):
    from atom_openmm.rbfe_structprep import rbfe_structprep
    from openmm import unit
    from openmm import app
    from openmm.app.amberprmtopfile import AmberPrmtopFile
    from openmm import XmlSerializer
    import numpy as np

    run_dir = os.path.join(tmp_path, "QB_A08_A07")
    shutil.copytree(os.path.join(curr_dir, "QB_A08_A07"), run_dir)

    box_vectors = np.array([[80.1356, 0, 0], [0, 86.0443, 0], [0, 0, 81.4926]])
    a = unit.Quantity(box_vectors[0] * unit.angstrom)
    b = unit.Quantity(box_vectors[1] * unit.angstrom)
    c = unit.Quantity(box_vectors[2] * unit.angstrom)
    box_vectors = (a, b, c)
    omm_structure = AmberPrmtopFile(
        os.path.join(run_dir, "QB_A08_A07.prmtop"), periodicBoxVectors=box_vectors
    )
    system = omm_structure.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=9 * unit.angstrom,
        constraints=app.HBonds,
        hydrogenMass=1.5 * unit.amu,
        rigidWater=True,
        switchDistance=7 * unit.angstrom,
    )
    with open(os.path.join(run_dir, "QB_A08_A07_sys.xml"), "w") as f:
        f.write(XmlSerializer.serialize(system))

    keywords = None
    with open("QB_A08_A07_asyncre.yaml") as f:
        keywords = yaml.safe_load(f)
        keywords['WORKDIR'] = run_dir
    rbfe_structprep(config_file = None, options = keywords)

def _test_rbfe_production(tmp_path):
    from atom_openmm.rbfe_production import rbfe_production

    run_dir = os.path.join(tmp_path, "QB_A08_A07_equil")
    shutil.copytree(os.path.join(curr_dir, "QB_A08_A07_equil"), run_dir)

    rbfe_production(os.path.join(run_dir, "QB_A08_A07_asyncre.yaml"))

def _test_rbfe_production_from_keywords(tmp_path):
    from atom_openmm.rbfe_production import rbfe_production

    run_dir = os.path.join(tmp_path, "QB_A08_A07_equil")
    shutil.copytree(os.path.join(curr_dir, "QB_A08_A07_equil"), run_dir)
    keywords = None
    with open("QB_A08_A07_asyncre.yaml") as f:
        keywords = yaml.safe_load(f)
        keywords['WORKDIR'] = run_dir
    rbfe_production(config_file = None, options = keywords)

def _test_make_atm_system_from_amber(tmp_path):
    from atom_openmm.make_atm_system_from_amber import make_system

    refxml = os.path.join(curr_dir, "QB_A08_A07_equil", "QB_A08_A07_sys.xml")
    outxml = os.path.join(tmp_path, "QB_A08_A07_sys.xml")
    prmtopfile = os.path.join(curr_dir, "QB_A08_A07", "QB_A08_A07.prmtop")
    crdfile = os.path.join(curr_dir, "QB_A08_A07", "QB_A08_A07.inpcrd")
    pdboutfile = os.path.join(tmp_path, "QB_A08_A07_sys.pdb")
    make_system(
        prmtopfile=prmtopfile,
        crdfile=crdfile,
        xmloutfile=outxml,
        pdboutfile=pdboutfile,
        hmass=1.5,
        switchDistance=0.7,
        nonbondedCutoff=0.9,
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    with open(refxml, "r") as f:
        ref_lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"

    os.remove(outxml)
    os.system(
        f"make_atm_system_from_amber --AmberPrmtopinFile {prmtopfile} --AmberInpcrdinFile {crdfile} --systemXMLoutFile {outxml} --systemPDBoutFile {pdboutfile} --hmass 1.5 --switchDistance 0.7 --nonbondedCutoff 0.9"
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"


def _test_make_atm_system_from_pdb(tmp_path):
    from atom_openmm.make_atm_system_from_pdb import make_system

    refxml = os.path.join(curr_dir, "3ptb", "3ptb_sys.xml")
    outxml = os.path.join(tmp_path, "3ptb_sys.xml")
    pdb = os.path.join(curr_dir, "3ptb", "3ptb.pdb")
    sdffile = os.path.join(curr_dir, "3ptb", "BEN_ideal.sdf")
    pdboutfile = os.path.join(tmp_path, "3ptb_sys.pdb")
    make_system(
        systempdbfile=pdb,
        ligandsdffile=sdffile,
        lig1resid=1,
        xmloutfile=outxml,
        pdboutfile=pdboutfile,
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    with open(refxml, "r") as f:
        ref_lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"

    os.remove(outxml)
    os.system(
        f"make_atm_system_from_pdb --systemPDBinFile {pdb} --ligandsSDFFile {sdffile} --LIG1resid 1 --systemXMLoutFile {outxml} --systemPDBoutFile {pdboutfile}"
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"


def _test_make_atm_system_from_rcpt_lig(tmp_path):
    from atom_openmm.make_atm_system_from_rcpt_lig import make_system

    refxml = os.path.join(curr_dir, "3ptb", "3ptb_sys_2.xml")
    outxml = os.path.join(tmp_path, "3ptb_sys_2.xml")
    pdb = os.path.join(curr_dir, "3ptb", "3ptb.pdb")
    sdffile = os.path.join(curr_dir, "3ptb", "BEN_ideal.sdf")
    pdboutfile = os.path.join(tmp_path, "3ptb_sys_2.pdb")
    make_system(
        receptorfile=pdb,
        lig1file=sdffile,
        displacement=[22, 22, 22],
        xmloutfile=outxml,
        pdboutfile=pdboutfile,
        ionicstrength=0,
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    with open(refxml, "r") as f:
        ref_lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"

    os.remove(outxml)
    os.system(
        f"make_atm_system_from_rcpt_lig --receptorinFile {pdb} --LIG1inFile {sdffile} --displacement '22.0 22.0 22.0' --systemXMLoutFile {outxml} --systemPDBoutFile {pdboutfile} --ionicStrength 0.0 "
    )
    with open(outxml, "r") as f:
        lines = f.readlines()[2:]  # Skipping OpenMM version header
    assert lines == ref_lines, f"Failed comparison of XML files: {outxml} != {refxml}"


def _test_input_parser():
    from atom_openmm.utils.config import parse_config

    config_yaml = parse_config(os.path.join(curr_dir, "QB_A08_A07_asyncre.yaml"))
    config_json = parse_config(os.path.join(curr_dir, "QB_A08_A07_asyncre.json"))
    config_cntl = parse_config(os.path.join(curr_dir, "QB_A08_A07_asyncre.cntl"))

    assert config_yaml == config_json
    assert config_yaml == config_cntl

    for key in config_yaml.keys():
        assert config_yaml[key] == config_json[key]
        assert config_yaml[key] == config_cntl[key]
        assert type(config_yaml[key]) is type(
            config_json[key]
        ), f"{key} {type(config_yaml[key])} {type(config_json[key])}"
        assert type(config_yaml[key]) is type(
            config_cntl[key]
        ), f"{key} {type(config_yaml[key])} {type(config_cntl[key])}"


def _test_sync_production_ats(tmp_path):
    from atom_openmm.rbfe_production import rbfe_production

    run_dir = os.path.join(tmp_path, "QB_A08_A07_equil_sync_ats")
    shutil.copytree(os.path.join(curr_dir, "QB_A08_A07_equil_sync_ats"), run_dir)

    rbfe_production(os.path.join(run_dir, "QB_A08_A07.yaml"))

    for i in range(4):
        assert os.path.exists(os.path.join(run_dir, f"r{i}", "QB_A08_A07.xtc"))


def _test_sync_production_incremental(tmp_path):
    from atom_openmm.rbfe_production import rbfe_production
    from openmm.app.internal.xtc_utils import read_xtc
    import yaml

    run_dir = os.path.join(tmp_path, "QB_A08_A07_equil_sync")
    shutil.copytree(os.path.join(curr_dir, "QB_A08_A07_equil_sync"), run_dir)

    configfile = os.path.join(run_dir, "QB_A08_A07.yaml")
    with open(configfile, "r") as f:
        config = yaml.safe_load(f)
    config["MAX_SAMPLES"] = "+1"
    with open(configfile, "w") as f:
        yaml.dump(config, f)

    rbfe_production(configfile)
    for i in range(4):
        xtcf = os.path.join(run_dir, f"r{i}", "QB_A08_A07.xtc")
        coords_read, box_read, time, step = read_xtc(xtcf.encode("utf-8"))
        assert coords_read.shape == (46380, 3, 1)
        assert box_read.shape == (3, 3, 1)
        assert len(time) == 1
        assert len(step) == 1

    startsampl_file = os.path.join(run_dir, "starting_sample")
    with open(startsampl_file, "r") as f:
        starting_sample = int(f.read().strip())
        assert starting_sample == 0
    prog_file = os.path.join(run_dir, "progress")
    with open(prog_file, "r") as f:
        progress = float(f.read().strip())
        assert progress == 0.0

    # Test restarting
    config["MAX_SAMPLES"] = "+2"
    with open(configfile, "w") as f:
        yaml.dump(config, f)
    rbfe_production(configfile)
    for i in range(4):
        xtcf = os.path.join(run_dir, f"r{i}", "QB_A08_A07.xtc")
        coords_read, box_read, time, step = read_xtc(xtcf.encode("utf-8"))
        assert coords_read.shape == (46380, 3, 2)
        assert box_read.shape == (3, 3, 2)
        assert len(time) == 2
        assert len(step) == 2

    with open(startsampl_file, "r") as f:
        starting_sample = int(f.read().strip())
        assert (
            starting_sample == 0
        )  # The starting_sample file existed so was not overwritten
    with open(prog_file, "r") as f:
        progress = float(f.read().strip())
        assert progress == 0.5

    # Add two more samples
    os.remove(startsampl_file)
    os.remove(prog_file)
    config["MAX_SAMPLES"] = "+2"
    with open(configfile, "w") as f:
        yaml.dump(config, f)
    rbfe_production(configfile)
    for i in range(4):
        xtcf = os.path.join(run_dir, f"r{i}", "QB_A08_A07.xtc")
        coords_read, box_read, time, step = read_xtc(xtcf.encode("utf-8"))
        assert coords_read.shape == (46380, 3, 4)
        assert box_read.shape == (3, 3, 4)
        assert len(time) == 4
        assert len(step) == 4

    with open(startsampl_file, "r") as f:
        starting_sample = int(f.read().strip())
        assert (
            starting_sample == 2
        )  # We removed the starting_sample file so it was written from the checkpoint
    with open(prog_file, "r") as f:
        progress = float(f.read().strip())
        assert progress == 0.5


def _test_xtc_writing_every_n_cycles(tmp_path):
    from atom_openmm.rbfe_production import rbfe_production
    from openmm.app.internal.xtc_utils import read_xtc
    import yaml

    run_dir = os.path.join(tmp_path, "QB_A08_A07_equil_sync")
    shutil.copytree(os.path.join(curr_dir, "QB_A08_A07_equil_sync"), run_dir)

    configfile = os.path.join(run_dir, "QB_A08_A07.yaml")
    with open(configfile, "r") as f:
        config = yaml.safe_load(f)
    config["MAX_SAMPLES"] = "3"
    config["PRODUCTION_STEPS"] = "10"
    config["TRJ_FREQUENCY"] = "20"
    with open(configfile, "w") as f:
        yaml.dump(config, f)

    rbfe_production(configfile)
    for i in range(4):
        xtcf = os.path.join(run_dir, f"r{i}", "QB_A08_A07.xtc")
        coords_read, box_read, time, step = read_xtc(xtcf.encode("utf-8"))
        assert coords_read.shape == (46380, 3, 1)
        assert box_read.shape == (3, 3, 1)
        assert len(time) == 1
        assert len(step) == 1


def _test_multi_var_atoms_atmforce(tmp_path):
    import copy
    import logging
    import numpy as np
    import openmm as mm
    from atom_openmm.ommsystem import OMMSystemRBFE

    def _assert_fails(fn):
        try:
            fn()
        except AssertionError:
            return
        raise AssertionError("Expected AssertionError for inconsistent atom lists")

    ats_dir = os.path.join(curr_dir, "QB_A08_A07_equil_sync_ats")
    with open(os.path.join(ats_dir, "QB_A08_A07.yaml"), "r") as f:
        base_kw = yaml.safe_load(f)
    pdbfile = os.path.join(ats_dir, "QB_A08_A07.pdb")
    sysfile = os.path.join(ats_dir, "QB_A08_A07_sys.xml")
    logger = logging.getLogger("test_multi_var")

    # 1. Particle-based displacements with 3 variable regions
    kw_p3 = copy.deepcopy(base_kw)
    kw_p3["LIGAND1_VAR_ATOMS"] = [[4671, 4673], [4682, 4683], [4684, 4685, 4686]]
    kw_p3["LIGAND2_VAR_ATOMS"] = [[4715, 4716, 4717], [4718, 4719, 4720], [4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728]]
    kw_p3["LIGAND1_ATTACH_ATOM"] = [4669, 4668, 4667]
    kw_p3["LIGAND2_ATTACH_ATOM"] = [4704, 4703, 4702]
    kw_p3.pop("LIGAND1_COMMON_ATOMS", None)
    kw_p3.pop("LIGAND2_COMMON_ATOMS", None)
    syst_p3 = OMMSystemRBFE("QB_A08_A07", kw_p3, pdbfile, sysfile, logger)
    syst_p3.create_system()

    # Inconsistent number of variable atom lists must fail
    kw_bad = copy.deepcopy(kw_p3)
    kw_bad["LIGAND2_VAR_ATOMS"] = [[4715, 4716], [4718, 4719]]
    _assert_fails(lambda: OMMSystemRBFE("QB_A08_A07", kw_bad, pdbfile, sysfile, logger).create_system())

    # 2. Frame-based displacements with 1, 2, 3, 4 lists and empty sublist
    if hasattr(mm, "ParticleFrameDisplacement"):
        base_frame = copy.deepcopy(base_kw)
        base_frame.pop("LIGAND1_ATTACH_ATOM", None)
        base_frame.pop("LIGAND2_ATTACH_ATOM", None)
        base_frame.pop("LIGAND1_COMMON_ATOMS", None)
        base_frame.pop("LIGAND2_COMMON_ATOMS", None)
        base_frame["VAR_ATOMS_MASS_SCALE"] = 2.0

        # 3 variable regions
        kw_f3 = copy.deepcopy(base_frame)
        kw_f3["LIGAND1_VAR_ATOMS"] = [[4671, 4673], [4682, 4683], [4684, 4685, 4686]]
        kw_f3["LIGAND2_VAR_ATOMS"] = [[4715, 4716, 4717], [4718, 4719], [4720, 4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728]]
        kw_f3["ALIGN_LIGAND1_REF_ATOMS"] = [[14, 11, 15], [17, 14, 18], [10, 9, 11]]
        kw_f3["ALIGN_LIGAND2_REF_ATOMS"] = [[14, 11, 15], [17, 14, 18], [10, 9, 11]]
        syst_f3 = OMMSystemRBFE("QB_A08_A07", kw_f3, pdbfile, sysfile, logger)
        syst_f3.create_system()
        ctx = mm.Context(syst_f3.system, syst_f3.integrator, mm.Platform.getPlatformByName("Reference"))
        ctx.setPositions(syst_f3.positions)
        u1, u0, ebias = syst_f3.atmforce.getPerturbationEnergy(ctx)
        assert np.isfinite(u1.value_in_unit(mm.unit.kilojoule_per_mole))
        assert np.isfinite(u0.value_in_unit(mm.unit.kilojoule_per_mole))
        del ctx

        # 4 variable regions including an empty sublist
        kw_f4 = copy.deepcopy(base_frame)
        kw_f4["LIGAND1_VAR_ATOMS"] = [[], [4671, 4673], [4682, 4683], [4684, 4685, 4686]]
        kw_f4["LIGAND2_VAR_ATOMS"] = [[4715], [4716, 4717], [4718, 4719], [4720, 4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728]]
        kw_f4["ALIGN_LIGAND1_REF_ATOMS"] = [[14, 11, 15], [17, 14, 18], [10, 9, 11], [5, 4, 6]]
        kw_f4["ALIGN_LIGAND2_REF_ATOMS"] = [[14, 11, 15], [17, 14, 18], [10, 9, 11], [5, 4, 6]]
        syst_f4 = OMMSystemRBFE("QB_A08_A07", kw_f4, pdbfile, sysfile, logger)
        syst_f4.create_system()

        # Inconsistent REF_ATOMS count (2 vs 3) must fail
        kw_bad_ref = copy.deepcopy(kw_f3)
        kw_bad_ref["ALIGN_LIGAND1_REF_ATOMS"] = [[14, 11, 15], [17, 14, 18]]
        kw_bad_ref["ALIGN_LIGAND2_REF_ATOMS"] = [[14, 11, 15], [17, 14, 18]]
        _assert_fails(lambda: OMMSystemRBFE("QB_A08_A07", kw_bad_ref, pdbfile, sysfile, logger).create_system())
