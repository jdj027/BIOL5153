#! /usr/bin/env python

# Monte Carlo for Estimation of Endogenous Fluorescence Signaling in Skin
# Jake Jones
# April/May 2019

import random
import math

# Literature Cited:
# Meglinski and Matcher, 2002, Physiol. Meas., 23, 741
# Rehman et al., 2016, Biomed. Optics Exp., 8(3), 1488
# Pena et al., 2005, Optics Exp., 13(16), 6268-6274
# Chen et al., 2006, Scanning, 28, 319-326

#-----------------------------------------------------------------------------#
# Skin Layer Thickness measurements:                                          |
#-----------------------------------------------------------------------------|
# Stratum corneum = 20 um                                                     |
# Epidermis = 80 um                                                           |
# Papillary Dermis = 150 um                                                   |
# Upper blood net dermis = 100 um                                             |
# Reticular dermis = 1500 um                                                  |
# Deep blood net dermis = 100 um                                              |
# Subcutaneous fat = 2000 um                                                  |
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Skin layer anisotropy measurements:                                         |
#-----------------------------------------------------------------------------|
# Stratum corneum = 0.86                                                      |
# Epidermis = 0.8                                                             |
# Papillary Dermis = 0.9                                                      |
# Upper blood net dermis = 0.95                                               |
# Reticular dermis = 0.8                                                      |
# Deep blood net dermis = 0.95                                                |
# Subcutaneous fat = 0.75                                                     |
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Skin layer u_s measurements (1/mm):                                         |
#-----------------------------------------------------------------------------|
# Stratum corneum = 100                                                       |
# Epidermis = 45                                                              |
# Papillary Dermis = 30                                                       |
# Upper blood net dermis = 35                                                 |
# Reticular dermis = 25                                                       |
# Deep blood net dermis = 30                                                  |
# Subcutaneous fat = 5                                                        |
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Skin layer blood percentage:                                                |
#-----------------------------------------------------------------------------|
# Stratum corneum = 0                                                         |
# Epidermis = 0                                                               |
# Papillary Dermis = 0.04                                                     |
# Upper blood net dermis = 0.3                                                |
# Reticular dermis = 0.04                                                     |
# Deep blood net dermis = 0.1                                                 |
# Subcutaneous fat = 0.05                                                     |
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Skin layer water percentage:                                                |
#-----------------------------------------------------------------------------|
# Stratum corneum = 0.05                                                      |
# Epidermis = 0.2                                                             |
# Papillary Dermis = 0.5                                                      |
# Upper blood net dermis = 0.6                                                |
# Reticular dermis = 0.7                                                      |
# Deep blood net dermis = 0.7                                                 |
# Subcutaneous fat = 0.7                                                      |
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Skin layer cellular density percentage:                                     |
#-----------------------------------------------------------------------------|
# Stratum corneum = 0                                                         |
# Epidermis = 0.975                                                           |
# Papillary Dermis = 0.9                                                      |
# Upper blood net dermis = 0.6                                                |
# Reticular dermis = 0.25                                                     |
# Deep blood net dermis = 0.4                                                 |
# Subcutaneous fat = 0.15                                                     |
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Quantum yield of endogenous fluorophores:                                   |
#-----------------------------------------------------------------------------|
# Free NADH = 0.02                                                            |
# Bound NADH = 0.1                                                            |
# Flavins = 0.25                                                              |
# Keratin = 0.65                                                              |
#-----------------------------------------------------------------------------#
# Assumptions:
    # White light source
    # Minimal melanin in sample
    # Collagen AF, crosslink AF are grouped into debris/connective tissue count


u_s = 90.0              # scattering coefficient (1/cm)
u_a = 10.0              # absorption coefficient (1/cm)
g = 0.75                # anisotropy factor (g=1 would be isotropic)
photoncount = 50000     # Number of photons you want to run through
n_in = 1.4              # Refractive index of bulk skin
n_out = 1               # Refractive index of air

# Thickness of individual skin layers in cm and relative depth
sc_thickness = 0.002
sc_depth = sc_thickness

epi_thickness = 0.008
epi_depth = sc_depth + epi_thickness

papderm_thickness = 0.015
papderm_depth = epi_depth + papderm_thickness

ubnderm_thickness = 0.01
ubnderm_depth = papderm_depth + ubnderm_thickness

retderm_thickness = 0.15
retderm_depth = ubnderm_depth + retderm_thickness

dbnderm_thickness = 0.01
dbnderm_depth = retderm_depth + dbnderm_thickness

subcut_thickness = 0.2
subcut_depth = dbnderm_depth + subcut_thickness

# Total thickness of skin "slab"
thickness = sc_thickness + epi_thickness + papderm_thickness + ubnderm_thickness + retderm_thickness + dbnderm_thickness + subcut_thickness

# Absorptive species percent concentrations in each layer per skin slab
sc_blood = 0.0
sc_water = 0.05
sc_abs = sc_blood + sc_water

epi_blood = 0.0
epi_water = 0.2
epi_abs = epi_blood + epi_water

papderm_blood = 0.04
papderm_water = 0.5
papderm_abs = papderm_blood + papderm_water

ubnderm_blood = 0.3
ubnderm_water = 0.6
ubnderm_abs = ubnderm_blood + ubnderm_water

retderm_blood = 0.04
retderm_water = 0.7
retderm_abs = retderm_blood + retderm_water

dbnderm_blood = 0.1
dbnderm_water = 0.7
dbnderm_abs = dbnderm_blood + dbnderm_water

subcut_blood = 0.05
subcut_water = 0.7
subcut_abs = subcut_blood + subcut_water

# Fluorescence Quantum Yields of Endogenous Fluorophores
qy_nadh = 0.12
qy_fad = 0.25
qy_keratn = 0.65

def speckle_ref(photoncount):
    speckle_effect = ((n_in - n_out) ** 2)/((n_in + n_out) ** 2)
    reflected_photons = int(round(speckle_effect * photoncount))
    incident_photons = photoncount - reflected_photons

    return incident_photons, reflected_photons

def photon_travel(u_a,u_s,g):
    # Initialize photon "weight": Weight is a simplied quantum mechanics consideration that can be used to determine when a photon "dies" based on its energy
    photon_weight = 1

    # Inital interaction coefficient (1/cm)
    u_t = u_s + u_a

    # Initialize photon position at surface of skin (Assumption: Only care about objective powers, i.e. the photons that make it to the sample)
    posx = 0
    posy = 0
    posz = 0

    # Initialize counting variables
    absorbed_chromo = 0
    transmitted_count = 0
    reflected_count = 0
    nadh_count = 0
    nadh_abs_count = 0
    fad_count = 0
    fad_abs_count = 0
    keratin_count = 0
    keratin_abs_count = 0
    debris_count = 0

    # Initialize directional cosines so that photon is entering normal to the tissue surface (Assumption: light source = laser)
    cox = 0
    coy = 0
    coz = 1

    while photon_weight > 0:
        random_step = (float(random.randint(1,1000))) / 1000
        step_size = (-1 * math.log(random_step)) / (float(u_a) + float(u_s))

        # Photon moves on to new coordinates reached by random step and directional cosines
        posx = posx + (cox*step_size)
        posy = posy + (coy*step_size)
        posz = posz + (coz*step_size)

        # Check if the photon is reflected, transmitted, or absorbed.
        # This is where we will need to check layers. The reflected and transmitted
        # statements are simple boundary condition based on the dimensions of the slab.
        # The absorbed section is the most complex and is dependent on skin layer.
        # Z-position decides the layer, then updates environmental variables and
        # establishes a blood/water/cell density before deciding what the photon
        # packet is absorbed by and if it causes a fluorescence emission.
        if (posz < 0) and (photon_weight > 0):
            reflected_count = reflected_count + photon_weight
            photon_weight = 0
            break
        elif (posz > thickness) and (photon_weight > 0):
            transmitted_count = transmitted_count + photon_weight
            photon_weight = 0
            break
        else:
            # Determine which layer of skin the photon interacts with
            # Once layer is determined, update: g, us, ut, chormophore %, cell %

            # Stratum Corneum
            if (posz >= 0) and (posz <= sc_depth) and (photon_weight > 0):
                g = 0.86
                u_s = 10.0
                u_t = u_a + u_s
                chrom_percent = sc_abs
                cell_percent = 0.0
                ###############################################################
                # need a unique check for keratin af
                ###############################################################
                #print('SC')
            # Epidermis
            elif (posz > sc_depth) and (posz <= epi_depth) and (photon_weight > 0):
                g = 0.8
                u_s = 4.5
                u_t = u_a + u_s
                chrom_percent = epi_abs
                cell_percent = 0.975
                #print('Epidermis')
            # Papillary Dermis
            elif (posz > epi_depth) and (posz <= papderm_depth) and (photon_weight > 0):
                g = 0.9
                u_s = 3.0
                u_t = u_a + u_s
                chrom_percent = papderm_abs
                cell_percent = 0.9
                #print('Pap Derm')
            # Upper Blood Net Dermis
            elif (posz > papderm_depth) and (posz <= ubnderm_depth) and (photon_weight > 0):
                g = 0.95
                u_s = 3.5
                u_t = u_a + u_s
                chrom_percent = ubnderm_abs
                cell_percent = 0.6
                #print('UBN Derm')
            # Reticular Dermis
            elif (posz > ubnderm_depth) and (posz <= retderm_depth) and (photon_weight > 0):
                g = 0.8
                u_s = 2.5
                u_t = u_a + u_s
                chrom_percent = retderm_abs
                cell_percent = 0.25
                #print('Reticular Derm')
            # Deep Blood Net Dermis
            elif (posz > retderm_depth) and (posz <= dbnderm_depth) and (photon_weight > 0):
                g = 0.95
                u_s = 3.0
                u_t = u_a + u_s
                chrom_percent = dbnderm_abs
                cell_percent = 0.4
                #print('DBN Derm')
            # Subcutaneous Fat
            elif (posz > dbnderm_depth) and (posz <= subcut_depth) and (photon_weight > 0):
                g = 0.7
                u_s = 0.5
                u_t = u_a + u_s
                chrom_percent = subcut_abs
                cell_percent = 0.15
                #print('Subcut')

            # RNG if the photon hits a cell, chromophore or connective tissue/debris
            absorption_interaction = (float(random.randint(1,1000))) / 1000
            cell_impact = (float(random.randint(1,1000))) / 1000

            # If cell is avoided, check for blood/water chromatic absorption
            if (absorption_interaction <= chrom_percent) and (cell_impact > cell_percent) and (photon_weight > 0):
                # chromophore or water, this happens
                albedo = float(u_a) / float(u_t)
                weight_absorbed = float(photon_weight) * albedo
                photon_weight = photon_weight - weight_absorbed
                absorbed_chromo = absorbed_chromo + weight_absorbed

                # Photon packet survival roulette to mimic quantum "death"
                if photon_weight < 0.005:
                    challenge = 20
                    roulette = (float(random.randint(1,1000))) / 1000
                    if roulette <= (1/challenge):
                        photon_weight = challenge*photon_weight # Survival = amplified weight reduction
                    else:
                        photon_weight = 0 # Death = end of photon packet travel
                        break

                # If the photon survives the roulette, it now scatters
                # Generate random #s to form scattering angles (polar and azimuthal)
                random_polarangle = (float(random.randint(1,1000))) / 1000
                random_azimuthalangle = (float(random.randint(1,1000))) / 1000

                # Calculate the polar scattering angle
                if g==0:
                    polar = math.acos(2 * random_polarangle - 1)
                else:
                    polar = math.acos((1/(2*g))*(1+(g*g)-(((1-(g*g))/((1-g)+(2*g*random_polarangle))))))

                # Calculate the azimuthal scattering angle
                azimuth = 2*math.pi*random_azimuthalangle

                # Update the directional cosines to change the photon packet vector
                if math.fabs(coz) > 0.999:
                    cox = math.sin(polar)*math.cos(azimuth)
                    coy = math.sin(polar)*math.sin(azimuth)
                    coz = math.copysign(coz,0.0)*math.cos(polar)
                else:
                    cox = math.sin(polar)*((cox*coz*math.cos(azimuth)-coy*math.sin(azimuth))/math.sqrt(1-(coz*coz)))+cox*math.cos(polar);
                    coy = math.sin(polar)*((coy*coz*math.cos(azimuth)+cox*math.sin(azimuth))/math.sqrt(1-(coz*coz)))+coy*math.cos(polar);
                    coz = -math.sin(polar)*math.cos(azimuth)*math.sqrt(1-(coz*coz))+coz*math.cos(polar);

            # If a photon hits a cell, determine if fluorescence is emitted
            elif (absorption_interaction > chrom_percent) and (cell_impact <= cell_percent) and (posz > sc_depth):
                fluorophore_coinflip = (float(random.randint(1,1000))) / 1000
                fluorescence_generation = (float(random.randint(1,1000))) / 1000

                albedo = float(u_a) / float(u_t)
                weight_absorbed = float(photon_weight) * albedo
                photon_weight = photon_weight - weight_absorbed

                # Determine if a fluorophore absorbs a photon, or creates a new fluorescent photon packet
                if (fluorophore_coinflip < 0.5):
                    if fluorescence_generation <= qy_fad:
                        fad_count = fad_count + weight_absorbed
                    else:
                        fad_abs_count = fad_abs_count + weight_absorbed
                else:
                    if fluorescence_generation <= qy_nadh:
                        nadh_count = nadh_count + weight_absorbed
                    else:
                        nadh_abs_count = nadh_abs_count + weight_absorbed

                # Photon packet survival roulette to mimic quantum "death"
                if photon_weight < 0.005:
                    challenge = 20
                    roulette = (float(random.randint(1,1000))) / 1000
                    if roulette <= (1/challenge):
                        photon_weight = challenge*photon_weight # Survival = amplified weight reduction
                    else:
                        photon_weight = 0 # Death = end of photon packet travel
                        break

                # If the photon survives the roulette, it now scatters
                # Generate random #s to form scattering angles (polar and azimuthal)
                random_polarangle = (float(random.randint(1,1000))) / 1000
                random_azimuthalangle = (float(random.randint(1,1000))) / 1000

                # Calculate the polar scattering angle
                if g==0:
                    polar = math.acos(2 * random_polarangle - 1)
                else:
                    polar = math.acos((1/(2*g))*(1+(g*g)-(((1-(g*g))/((1-g)+(2*g*random_polarangle))))))

                # Calculate the azimuthal scattering angle
                azimuth = 2*math.pi*random_azimuthalangle

                # Update the directional cosines to change the photon packet vector
                if math.fabs(coz) > 0.999:
                    cox = math.sin(polar)*math.cos(azimuth)
                    coy = math.sin(polar)*math.sin(azimuth)
                    coz = math.copysign(coz,0.0)*math.cos(polar)
                else:
                    cox = math.sin(polar)*((cox*coz*math.cos(azimuth)-coy*math.sin(azimuth))/math.sqrt(1-(coz*coz)))+cox*math.cos(polar);
                    coy = math.sin(polar)*((coy*coz*math.cos(azimuth)+cox*math.sin(azimuth))/math.sqrt(1-(coz*coz)))+coy*math.cos(polar);
                    coz = -math.sin(polar)*math.cos(azimuth)*math.sqrt(1-(coz*coz))+coz*math.cos(polar);

            # If the photon misses cells and does not interact with blood/water then it hits connective tissue/debris
            else:
                albedo = float(u_a) / float(u_t)
                weight_absorbed = float(photon_weight) * albedo
                photon_weight = photon_weight - weight_absorbed

                # Photon interacts with collagen/connective tissue or other debris
                debris_count = debris_count + photon_weight

                # Photon packet survival roulette to mimic quantum "death"
                if photon_weight < 0.005:
                    challenge = 20
                    roulette = (float(random.randint(1,1000))) / 1000
                    if roulette <= (1/challenge):
                        photon_weight = challenge*photon_weight # Survival = amplified weight reduction
                    else:
                        photon_weight = 0 # Death = end of photon packet travel
                        break

                # If the photon survives the roulette, it now scatters
                # Generate random #s to form scattering angles (polar and azimuthal)
                random_polarangle = (float(random.randint(1,1000))) / 1000
                random_azimuthalangle = (float(random.randint(1,1000))) / 1000

                # Calculate the polar scattering angle
                if g==0:
                    polar = math.acos(2 * random_polarangle - 1)
                else:
                    polar = math.acos((1/(2*g))*(1+(g*g)-(((1-(g*g))/((1-g)+(2*g*random_polarangle))))))

                # Calculate the azimuthal scattering angle
                azimuth = 2*math.pi*random_azimuthalangle

                # Update the directional cosines to change the photon packet vector
                if math.fabs(coz) > 0.999:
                    cox = math.sin(polar)*math.cos(azimuth)
                    coy = math.sin(polar)*math.sin(azimuth)
                    coz = math.copysign(coz,0.0)*math.cos(polar)
                else:
                    cox = math.sin(polar)*((cox*coz*math.cos(azimuth)-coy*math.sin(azimuth))/math.sqrt(1-(coz*coz)))+cox*math.cos(polar);
                    coy = math.sin(polar)*((coy*coz*math.cos(azimuth)+cox*math.sin(azimuth))/math.sqrt(1-(coz*coz)))+coy*math.cos(polar);
                    coz = -math.sin(polar)*math.cos(azimuth)*math.sqrt(1-(coz*coz))+coz*math.cos(polar);

    return absorbed_chromo, transmitted_count, reflected_count, fad_count, nadh_count, fad_abs_count, nadh_abs_count, debris_count

def main():
    (incident_photons, speckle_reflected_photons) = speckle_ref(photoncount)
    # Initialize counting variables
    abs_count = 0
    trs_count = 0
    ref_count = 0
    nf_count = 0
    na_count = 0
    ff_count = 0
    fa_count = 0
    #keratin_count = 0
    #keratin_abs_count = 0
    ct_count = 0
    for photon in range(0,incident_photons):
        (absorbed_photons, transmitted_photons, reflected_photons, fad_photons, nadh_photons, fad_abs_photons, nadh_abs_photons, debris_photons) = photon_travel(u_a,u_s,g)
        abs_count = abs_count + absorbed_photons
        trs_count = trs_count + transmitted_photons
        ref_count = ref_count + reflected_photons
        ff_count = ff_count + fad_photons
        fa_count = fa_count + fad_abs_photons
        nf_count = nf_count + nadh_photons
        na_count = na_count + nadh_abs_photons
        ct_count = ct_count + debris_photons
    total_weight = abs_count + trs_count + ref_count + ff_count + fa_count + nf_count + na_count + ct_count + speckle_reflected_photons
    transmission_correction = incident_photons - total_weight
    print('Chromatically Absorbed count: ' + str(abs_count))
    print('Transmitted count: ' + str(trs_count + transmission_correction))
    print('Reflected count: ' + str(ref_count + speckle_reflected_photons))
    print('FAD Fluorescence: ' + str(ff_count))
    print('FAD Absorption: ' + str(fa_count))
    print('NADH Fluorescence: ' + str(nf_count))
    print('NADH Absorption: ' + str(na_count))
    print('Connective Tissue Absorption: ' + str(ct_count))
    print('Total Weight: ' + str(total_weight + transmission_correction + speckle_reflected_photons))
# Get the argument before calling main
# args = get_args()

# Execute the program by calling main
if __name__=="__main__":
    main()
