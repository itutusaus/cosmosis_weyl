from builtins import range
from scipy.interpolate import interp1d
from cosmosis.datablock import names, option_section
import sys


def setup(options):
    perbin = options.get_bool(option_section, "perbin", True)
    auto_only = options.get_bool(option_section, "auto_only", False)
    apply_to_cl = options.get_bool(option_section, "apply_to_cl", True)
    if apply_to_cl:
        print("Applying b bin biases to C_ell values")
    else:
        print("Applying b bin biases to xi values")
    return perbin, auto_only, apply_to_cl

def bins_count_cl(block):
    n_z_bins_shear = 0
    n_z_bins_pos = 0
    
    if block.has_section('galaxy_shear_cl'):
        if "galaxy_shear_xi" in block:
            raise ValueError("""You've asked to apply binwise bias to the C(l)s, but we found xis
                in the block, so this almost definitely unintended""")
        n_z_bins_shear = block["galaxy_shear_cl", "nbin_b"]
        n_z_bins_pos = block["galaxy_shear_cl", "nbin_a"]

    if block.has_section('galaxy_cl'):
        if "galaxy_xi" in block:
            raise ValueError("""You've asked to apply binwise bias to the C(l)s, but we found xis
                in the block, so this almost definitely unintended""")
        n_z_bins_pos = block["galaxy_cl", "nbin"]

    if n_z_bins_pos==0:
        raise ValueError("Used bin bias module with apply_to_cl=T but did not find any C_ell density tracers")

    return n_z_bins_pos, n_z_bins_shear

def bins_count_xi(block):
    n_z_bins_shear = 0
    n_z_bins_pos = 0
    
    if block.has_section('galaxy_shear_xi'):
        n_z_bins_shear = block["galaxy_shear_xi", "nbin_b"]
        n_z_bins_pos = block["galaxy_shear_xi", "nbin_a"]

    if block.has_section('galaxy_xi'):
        n_z_bins_pos = block["galaxy_xi", "nbin_a"]

    if n_z_bins_pos==0:
        raise ValueError("Used bin bias module with apply_to_cl=F but did not find any xi density tracers")

    return n_z_bins_pos, n_z_bins_shear


def execute(block, options):
    perbin, auto_only, apply_to_cl = options

    if apply_to_cl:
        n_z_bins_pos, n_z_bins_shear = bins_count_cl(block)
    else:
        n_z_bins_pos, n_z_bins_shear = bins_count_xi(block)


    # We may be doing per-bin biases or a single global value
    if perbin:
        # per-bin - use b1,b2,b3, ...
        biases = [block["b_bin_bias", "bhat%d" % pos_bin]
                  for pos_bin in range(1, n_z_bins_pos + 1)]
    else:
        # all the same - just use b0
        biases = [block["b_bin_bias", "bhat0"] for pos_bin in range(n_z_bins_pos)]

    z=block[names.matter_power_nl,"z"]
    sig8_z=block[names.growth_parameters,"sigma_8"]
    sig8z_fun=interp1d(z,sig8_z)
    sig8z=[sig8z_fun(0.3),sig8z_fun(0.475),sig8z_fun(0.625),sig8z_fun(0.775),sig8z_fun(0.9),sig8z_fun(1.0)]
        
    for pos_bin1 in range(n_z_bins_pos):
        bias1 = biases[pos_bin1]
        
        if apply_to_cl:
            if block.has_section('galaxy_shear_cl'):
                for shear_bin in range(n_z_bins_shear):
                    name = "bin_{}_{}".format(pos_bin1 + 1, shear_bin + 1)
                    block["galaxy_shear_cl", name] *= bias1

            if block.has_section('galaxy_intrinsic_cl'):
                for shear_bin in range(n_z_bins_shear):
                    name = "bin_{}_{}".format(pos_bin1 + 1, shear_bin + 1)
                    block["galaxy_intrinsic_cl", name] *= bias1/sig8z[pos_bin1]
                    #Dividing by sigma8(zeff) in GR to keep the GR galaxy_intrinsic_cls
                
            if block.has_section('galaxy_cl'):
                for pos_bin2 in range(pos_bin1 + 1):
                    if auto_only:
                        if pos_bin2 != pos_bin1:
                            continue
                    bias2 = biases[pos_bin2]
                    name = "bin_{}_{}".format(pos_bin1 + 1, pos_bin2 + 1)
                    block['galaxy_cl', name] *= bias1 * bias2

        else:
            if block.has_section('galaxy_xi'):
                for pos_bin2 in range(pos_bin1 + 1):
                    if auto_only:
                        if pos_bin2 != pos_bin1:
                            continue
                    bias2 = biases[pos_bin2]
                    name = "bin_{}_{}".format(pos_bin1 + 1, pos_bin2 + 1)
                    block['galaxy_xi', name] *= bias1 * bias2

            if block.has_section('galaxy_shear_xi'):
                for shear_bin in range(n_z_bins_shear):
                    name = "bin_{}_{}".format(pos_bin1 + 1, shear_bin + 1)
                    block['galaxy_shear_xi', name] *= bias1


    return 0
