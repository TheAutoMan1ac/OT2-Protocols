from opentrons import protocol_api
import math

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

metadata = {
    'protocolName': 'Plasmid Purification with NucleoMag Plasmid Kit',
    'author': 'Markus Eder',
    'description': 'Protocol for purifying plasmids from bacterial cultures using magnetic beads on a 96-well plate'
}
requirements = {"robotType": "OT-2", "apiLevel": "2.19"}

def add_parameters(parameters: protocol_api.Parameters):
    parameters.add_int(
        variable_name="sample_count",
        display_name="Number of samples",
        description="How many samples should be purified.",
        default=1,
        minimum=1,
        maximum=24
    )
    parameters.add_bool(
        variable_name="debug",
        display_name="Debugging mode",
        description="Enable/Disable debugging mode",
        default=False
    )

def run(protocol: protocol_api.ProtocolContext):
    tips_300_per_sample = 8
    tips_1000_per_sample = 7
    
    tips_300_count = protocol.params.sample_count * tips_300_per_sample
    tips_1000_count = protocol.params.sample_count * tips_1000_per_sample
    
    tips_racks_300_count = math.ceil(tips_300_count / 96)
    tips_racks_1000_count = math.ceil(tips_1000_count / 96)
    
    tipracks_300_slots = [5, 8, 11]
    tipracks_1000_slots = [4, 7, 10]
    
    # tipracks
    tipracks_300 = [protocol.load_labware(load_name='opentrons_96_tiprack_300ul', location=slot) for slot in tipracks_300_slots[:tips_racks_300_count]]
    tipracks_1000 = [protocol.load_labware(load_name='opentrons_96_tiprack_1000ul', location=slot) for slot in tipracks_1000_slots[:tips_racks_1000_count]]
    
    # pipettes
    p300_single = protocol.load_instrument(instrument_name='p300_single', mount='right', tip_racks=tipracks_300)
    p1000_single = protocol.load_instrument(instrument_name='p1000_single', mount='left', tip_racks=tipracks_1000)
    
    #labware reagents
    tuberack_reagents = protocol.load_labware(load_name='opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical', location=tuberack_reagents_slot)
    tuberack_reagents_slot = 6
    
    reagent_A1_liquid = protocol.define_liquid(name="A1", description="Resuspension Buffer", display_color="#FF0000")
    reagent_A2_liquid = protocol.define_liquid(name="A2", description="Lysis Buffer", display_color="#FF0000")
    reagent_S3_liquid = protocol.define_liquid(name="S3", description="Neutralization Buffer", display_color="#FF0000")
    reagent_PAB_liquid = protocol.define_liquid(name="PAB", description="Binding Buffer", display_color="#FF0000")
    reagent_ERB_liquid = protocol.define_liquid(name="ERB", description="Detoxification Buffer", display_color="#FF0000")
    reagent_AQ_liquid = protocol.define_liquid(name="AQ", description="Wash Buffer", display_color="#FF0000")
    reagent_AE_liquid = protocol.define_liquid(name="AE", description="Elution Buffer", display_color="#FF0000")
    c_beads_liquid = protocol.define_liquid(name="C-Beads", description="NucleoMag Clearing Beads", display_color="#FF0000")
    m_beads_liquid = protocol.define_liquid(name="M-Beads", description="NucleoMag M-Beads", display_color="#FF0000")
    
    reagent_A1_slot = 'A1'
    reagent_A2_slot = 'A2'
    reagent_S3_slot = 'B1'
    reagent_PAB_slot = 'A3'
    reagent_ERB_slot = 'B3'
    reagent_AQ_slot = 'B4'
    reagent_AE_slot = 'B2'
    reagent_c_beads_slot = 'C1'
    reagent_m_beads_slot = 'C2'
    
    tuberack_reagents[reagent_A1_slot].load_liquid(liquid=reagent_A1_liquid, volume=math.ceil((90 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_A2_slot].load_liquid(liquid=reagent_A2_liquid, volume=math.ceil((120 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_S3_slot].load_liquid(liquid=reagent_S3_liquid, volume=math.ceil((120 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_PAB_slot].load_liquid(liquid=reagent_PAB_liquid, volume=math.ceil((390 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_ERB_slot].load_liquid(liquid=reagent_ERB_liquid, volume=math.ceil((1800 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_AQ_slot].load_liquid(liquid=reagent_AQ_liquid, volume=math.ceil((1800 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_AE_slot].load_liquid(liquid=reagent_AE_liquid, volume=math.ceil((100 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_c_beads_slot].load_liquid(liquid=c_beads_liquid, volume=35 * protocol.params.sample_count * 1.1)
    tuberack_reagents[reagent_m_beads_slot].load_liquid(liquid=m_beads_liquid, volume=20 * protocol.params.sample_count * 1.1)
    
    #labware samples
    eppiracks_count = math.ceil(protocol.params.sample_count / 12)
    eppiracks_samples_slots = [2, 3]
    eppiracks_samples = [protocol.load_labware(load_name='opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', location=slot) for slot in eppiracks_samples_slots[:eppiracks_count]]
    samples_unpurified_slots = ['A1','A2','A3','A4','A5','A6','B1','B2','B3','B4','B5','B6',
                                'A1','A2','A3','A4','A5','A6','B1','B2','B3','B4','B5','B6']
    sample_unpurified_eppirack_slots = [0,0,0,0,0,0,0,0,0,0,0,0,
                                        1,1,1,1,1,1,1,1,1,1,1,1]
    samples_purified_slots = ['C1','C2','C3','C4','C5','C6','D1','D2','D3','D4','D5','D6',
                              'C1','C2','C3','C4','C5','C6','D1','D2','D3','D4','D5','D6']
    sample_purified_eppirack_slots = [0,0,0,0,0,0,0,0,0,0,0,0,
                                      1,1,1,1,1,1,1,1,1,1,1,1]
    
    # modules
    # magdeck = Magnetic Module GEN1
    mag_module = protocol.load_module('magdeck', 1)
    mag_plate_96well = mag_module.load_labware('nest_96_wellplate_2ml_deep')
    mag_plate_process1_slots = ['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12',
                                'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12']
    mag_plate_process2_slots = ['C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12',
                                'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12']

    if protocol.params.debug:
        DELAY_LYSIS = 0
        DELAY_CLEARING_BEADS_INCUBATE = 0
        DELAY_CLEARING_BEADS_ENGAGE = 0
        DELAY_WASH_ENGAGE = 0
        
        DELAY_RESUSPEND_M_BEADS = 0
        DELAY_RESUSPEND_WASH = 0
        DELAY_RESUSPEND_DRY_BEADS = 0

        DELAY_DRY_BEADS = 0


    for i in range(protocol.params.sample_count):
        # Step 1: Add 90uL of reagent A1 and mix 
        p300_single.transfer(90, tuberack_reagents[reagent_A1_slot], eppiracks_samples[sample_unpurified_eppirack_slots[i]][samples_unpurified_slots[i]], mix_after=(MIX_TIMES_THOROUGH, 50))
        p300_single.transfer(90, eppiracks_samples[sample_unpurified_eppirack_slots[i]][i], mag_plate_96well[mag_plate_process1_slots[i]])

        # Step 2: Add 120uL of reagent A2 and wait for 2-3min
        p300_single.transfer(120, tuberack_reagents[reagent_A2_slot], mag_plate_96well[mag_plate_process1_slots[i]], mix_after=(MIX_TIMES_DEFAULT, 100))
        protocol.delay(minutes=DELAY_LYSIS)

        # Step 3: Add 120uL reagent S3 and mix
        p300_single.transfer(120, tuberack_reagents[reagent_S3_slot], mag_plate_96well[mag_plate_process1_slots[i]], mix_after=(MIX_TIMES_DEFAULT, 150))

        # Step 4: Add 35uL NucleoMag Clearing Beads and mix
        p300_single.transfer(35, tuberack_reagents[reagent_c_beads_slot], mag_plate_96well[mag_plate_process1_slots[i]], mix_before=(MIX_TIMES_DEFAULT, 20), mix_after=(MIX_TIMES_THOROUGH, 200))
        # not specified in protocol but feels right
        protocol.delay(minutes=DELAY_CLEARING_BEADS_INCUBATE)

        # Step 5: Magnetic separation
        mag_module.engage(height_from_base=ENGAGE_HEIGHT_DEFAULT)
        protocol.delay(minutes=DELAY_CLEARING_BEADS_ENGAGE)
        # Step 6: Transfer supernatant to new well
        p1000_single.transfer(365, mag_plate_96well[mag_plate_process1_slots[i]], mag_plate_96well[mag_plate_process2_slots[i]])
        mag_module.disengage()

        # Step 7: Add 20uL of NucleoMag M-Beads and 390uL of reagent PAB, mix, and wait
        p300_single.transfer(20, reagent_m_beads_slot, mag_plate_96well[mag_plate_process2_slots[i]], mix_before=(MIX_TIMES_DEFAULT, 10))
        p1000_single.transfer(390, tuberack_reagents[reagent_PAB_slot], mag_plate_96well[mag_plate_process2_slots[i]], mix_after=(MIX_TIMES_THOROUGH, 400))
        protocol.delay(minutes=DELAY_RESUSPEND_M_BEADS)

        # Step 8: Magnetic separation and remove supernatant
        mag_module.engage(height_from_base=ENGAGE_HEIGHT_DEFAULT)
        protocol.delay(minutes=DELAY_WASH_ENGAGE)
        p1000_single.pick_up_tip()
        p1000_single.aspirate(775, mag_plate_96well[mag_plate_process2_slots[i]])  # waste
        p1000_single.drop_tip()
        mag_module.disengage()

        # Steps 9-11: Wash with 900uL of ERB and AQ reagent, mix, remove supernatant and repeat
        # Wash 1 and 2 with ERB; Wash 3 and 4 with AQ 
        for i in range(4):
            p1000_single.pick_up_tip()
            if i == 0 or i == 1:
                p1000_single.transfer(900, tuberack_reagents[reagent_ERB_slot], mag_plate_96well[mag_plate_process2_slots[i]], mix_after=(MIX_TIMES_THOROUGH, 500), new_tip="never")
            elif i == 2 or i == 3:
                p1000_single.transfer(900, tuberack_reagents[reagent_AQ_slot], mag_plate_96well[mag_plate_process2_slots[i]], mix_after=(MIX_TIMES_THOROUGH, 500), new_tip="never")
                
            protocol.delay(minutes=DELAY_RESUSPEND_WASH)
            mag_module.engage(height_from_base=ENGAGE_HEIGHT_DEFAULT)
            protocol.delay(minutes=DELAY_WASH_ENGAGE)
            p1000_single.aspirate(900, mag_plate_96well[mag_plate_process2_slots[i]])  # waste
            p1000_single.drop_tip()
            mag_module.disengage()


        # Step 16: Let beads dry for 15min
        protocol.delay(minutes=DELAY_DRY_BEADS)
        
        # Step 17: Add 100uL of reagent AE, mix, and resuspend
        p300_single.transfer(100, tuberack_reagents[reagent_AE_slot], mag_plate_96well[mag_plate_process2_slots[i]], mix_after=(MIX_TIMES_THOROUGH, 50))

        # Step 18: Magnetic separation and transfer supernatant to eppendorf
        mag_module.engage(height_from_base=ENGAGE_HEIGHT_DEFAULT)
        protocol.delay(minutes=DELAY_RESUSPEND_DRY_BEADS)
        p300_single.transfer(100, mag_plate_96well[mag_plate_process2_slots[i]], eppiracks_samples[sample_purified_eppirack_slots[i]][samples_purified_slots[i]])  # final sample
        mag_module.disengage()