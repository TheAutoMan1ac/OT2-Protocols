from opentrons import protocol_api

# TODO?
# COULD ADD TIMER TO REGARD MIXING TIMES IN RESUSPENSION TIME

metadata = {
    'protocolName': 'Plasmid Purification with NucleoMag Plasmid Kit',
    'author': 'Markus Eder',
    'description': 'Protocol for purifying plasmids from bacterial cultures using magnetic beads on a 96-well plate'
}
requirements = {"robotType": "OT-2", "apiLevel": "2.19"}

def run(protocol: protocol_api.ProtocolContext):
    # modules
    # magdeck = Magnetic Module GEN1
    mag_module = protocol.load_module('magdeck', 1)
    
    # list of labware that meets the maximum volume requirements of the protocol (900uL)
    # nest_96_wellplate_2ml_deep
    # thermoscientificnunc_96_wellplate_1300ul
    # thermoscientificnunc_96_wellplate_2000ul
    # usascientific_96_wellplate_2.4ml_deep
    mag_plate = mag_module.load_labware('nest_96_wellplate_2ml_deep')

    # labware
    tiprack_300 = protocol.load_labware('opentrons_96_tiprack_300ul', 2)
    tiprack_1000 = protocol.load_labware('opentrons_96_tiprack_1000ul', 3)
    eppendorf_1p5 = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', 4)

    # pipettes
    # CHECK GENERATIONS AVAILABLE!
    p300_single = protocol.load_instrument('p300_single_gen2', mount='right', tip_racks=[tiprack_300])
    p1000_single = protocol.load_instrument('p1000_single_gen2', mount='left', tip_racks=[tiprack_1000])
    
    # reagents wells
    reagent_A1 = eppendorf_1p5['A1']
    reagent_A2 = eppendorf_1p5['A2']
    reagent_S3 = eppendorf_1p5['A3']
    clearing_beads = eppendorf_1p5['A4']
    nucleomag_m_beads = eppendorf_1p5['A5']
    reagent_PAB = eppendorf_1p5['A6']
    reagent_ERB = eppendorf_1p5['B1']
    reagent_AQ = eppendorf_1p5['B2']
    reagent_AE = eppendorf_1p5['B3']

    # processing and destination wells
    magplate_well_1 = mag_plate['A1']
    magplate_well_2 = mag_plate['A2']
    well_purified_plasmid = eppendorf_1p5['C1']

    DEBUG = False

    ENGAGE_HEIGHT_DEFAULT = 3

    # used for all mixing steps
    MIX_TIMES_DEFAULT = 5
    MIX_TIMES_THOROUGH = 10
    
    # minutes
    DELAY_LYSIS = 2 
    DELAY_CLEARING_BEADS_INCUBATE = 1
    DELAY_CLEARING_BEADS_ENGAGE = 1
    DELAY_WASH_ENGAGE = 1
    
    DELAY_RESUSPEND_M_BEADS = 5
    DELAY_RESUSPEND_WASH = 5
    DELAY_RESUSPEND_DRY_BEADS = 5

    DELAY_DRY_BEADS = 15

    if DEBUG:
        DELAY_LYSIS = 0
        DELAY_CLEARING_BEADS_INCUBATE = 0
        DELAY_CLEARING_BEADS_ENGAGE = 0
        DELAY_WASH_ENGAGE = 0
        
        DELAY_RESUSPEND_M_BEADS = 0
        DELAY_RESUSPEND_WASH = 0
        DELAY_RESUSPEND_DRY_BEADS = 0

        DELAY_DRY_BEADS = 0

    # Step 1: Add 90uL of reagent A1 and mix
    p300_single.transfer(90, reagent_A1, magplate_well_1, mix_after=(MIX_TIMES_THOROUGH, 50))

    # Step 2: Add 120uL of reagent A2 and wait for 2-3min
    p300_single.transfer(120, reagent_A2, magplate_well_1, mix_after=(MIX_TIMES_DEFAULT, 100))
    protocol.delay(minutes=DELAY_LYSIS)

    # Step 3: Add 120uL reagent S3 and mix
    p300_single.transfer(120, reagent_S3, magplate_well_1, mix_after=(MIX_TIMES_DEFAULT, 150))

    # Step 4: Add 35uL NucleoMag Clearing Beads and mix
    p300_single.transfer(35, clearing_beads, magplate_well_1, mix_before=(MIX_TIMES_DEFAULT, 20), mix_after=(MIX_TIMES_THOROUGH, 200))
    # not specified in protocol but feels right
    protocol.delay(minutes=DELAY_CLEARING_BEADS_INCUBATE)

    # Step 5: Magnetic separation
    mag_module.engage(height_from_base=ENGAGE_HEIGHT_DEFAULT)
    protocol.delay(minutes=DELAY_CLEARING_BEADS_ENGAGE)
    # Step 6: Transfer supernatant to new well
    p1000_single.transfer(365, magplate_well_1, magplate_well_2)
    mag_module.disengage()

    # Step 7: Add 20uL of NucleoMag M-Beads and 390uL of reagent PAB, mix, and wait
    p300_single.transfer(20, nucleomag_m_beads, magplate_well_2, mix_before=(MIX_TIMES_DEFAULT, 10))
    p1000_single.transfer(390, reagent_PAB, magplate_well_2, mix_after=(MIX_TIMES_THOROUGH, 400))
    protocol.delay(minutes=DELAY_RESUSPEND_M_BEADS)

    # Step 8: Magnetic separation and remove supernatant
    mag_module.engage(height_from_base=ENGAGE_HEIGHT_DEFAULT)
    protocol.delay(minutes=DELAY_WASH_ENGAGE)
    p1000_single.pick_up_tip()
    p1000_single.aspirate(775, magplate_well_2)  # waste
    p1000_single.drop_tip()
    mag_module.disengage()

    # Steps 9-11: Wash with 900uL of ERB and AQ reagent, mix, remove supernatant and repeat
    # Wash 1 and 2 with ERB; Wash 3 and 4 with AQ 
    for i in range(4):
        if i == 0 or i == 1:
            p1000_single.transfer(900, reagent_ERB, magplate_well_2, mix_after=(MIX_TIMES_THOROUGH, 500))
        elif i == 2 or i == 3:
            p1000_single.transfer(900, reagent_AQ, magplate_well_2, mix_after=(MIX_TIMES_THOROUGH, 500))
            
        protocol.delay(minutes=DELAY_RESUSPEND_WASH)
        mag_module.engage(height_from_base=ENGAGE_HEIGHT_DEFAULT)
        protocol.delay(minutes=DELAY_WASH_ENGAGE)
        p1000_single.pick_up_tip()
        p1000_single.aspirate(900, magplate_well_2)  # waste
        p1000_single.drop_tip()
        mag_module.disengage()


    # Step 16: Let beads dry for 15min
    protocol.delay(minutes=DELAY_DRY_BEADS)
    
    # Step 17: Add 100uL of reagent AE, mix, and resuspend
    p300_single.transfer(100, reagent_AE, magplate_well_2, mix_after=(MIX_TIMES_THOROUGH, 50))

    # Step 18: Magnetic separation and transfer supernatant to eppendorf
    mag_module.engage(height_from_base=ENGAGE_HEIGHT_DEFAULT)
    protocol.delay(minutes=DELAY_RESUSPEND_DRY_BEADS)
    p300_single.transfer(100, magplate_well_2, well_purified_plasmid)  # final sample
    mag_module.disengage()