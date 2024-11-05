from opentrons import protocol_api
import math
import time
import heapq

metadata = {
    'protocolName': 'Plasmid Purification with NucleoMag Plasmid Kit (Magnetic Separation)',
    'author': 'Markus Eder',
    'description': 'Protocol for purifying plasmids from pelleted bacterial cultures using NucleoMag Plasmid Kit with magnetic separation on a Brandt 96-well plate.'
}
requirements = {"robotType": "OT-2", "apiLevel": "2.20"}

def add_parameters(parameters: protocol_api.Parameters):
    parameters.add_bool(
        variable_name="debug",
        display_name="Debugging mode",
        description="Enable/Disable debugging mode",
        default=False
    )
    parameters.add_int(
        variable_name="sample_count",
        display_name="Number of samples",
        description="How many samples should be purified.",
        default=1,
        minimum=1,
        maximum=24
    )
    parameters.add_int(
        variable_name="engage_height",
        display_name="Magnetic module engage height",
        description="Height from origin magnetic module magnets should engage to.",
        default=3,
        minimum=-2,
        maximum=4,
        unit='mm'
    )
    
    
    #PARAMETERS MIXING
    parameters.add_int(
        variable_name="mix_times_default",
        display_name="Mix samples default times",
        description="Default number of mixing cycles.",
        default=5,
        minimum=1,
        maximum=50,
        unit='cycles'
    )
    parameters.add_int(
        variable_name="mix_times_thorough",
        display_name="Mix samples thorough times",
        description="Thorough number of mixing cycles.",
        default=20,
        minimum=1,
        maximum=50,
        unit='cycles'
    )
    parameters.add_int(
        variable_name="mix_times_beads",
        display_name="Mix bead stocks times",
        description="Number of mixing cycles for bead stocks.",
        default=10,
        minimum=1,
        maximum=50,
        unit='cycles'
    )
    parameters.add_int(
        variable_name="mix_times_resuspend_culture",
        display_name="Mix pellet resuspend",
        description="Number of mixing cycles for bacterial pellet resuspension.",
        default=30,
        minimum=1,
        maximum=50,
        unit='cycles'
    )
    
    #PARAMETERS DELAY TIMES
    parameters.add_int(
        variable_name="delay_lysis",
        display_name="Delay lysis",
        description="Delay protocol for lysis step",
        default=3,
        minimum=1,
        maximum=5,
        unit='min'
    )
    parameters.add_int(
        variable_name="delay_cbeads_incubate",
        display_name="Delay C-Beads incubate",
        description="Delay protocol for incubation of C-beads.",
        default=1,
        minimum=1,
        maximum=10,
        unit='min'
    )
    parameters.add_int(
        variable_name="delay_separate_default",
        display_name="Delay separate default",
        description="Delay protocol for separation of C/M-beads.",
        default=1,
        minimum=1,
        maximum=5,
        unit='min'
    )
    parameters.add_int(
        variable_name="delay_resuspend_thorough",
        display_name="Delay M-Beads resuspend",
        description="Delay protocol for resuspension of magnetic beads.",
        default=5,
        minimum=1,
        maximum=10,
        unit='min'
    )
    parameters.add_int(
        variable_name="delay_resuspend_wash",
        display_name="Delay M-Beads resuspend wash",
        description="Delay protocol for resuspension of magnetic beads during washing steps.",
        default=2,
        minimum=1,
        maximum=10,
        unit='min'
    )
    parameters.add_int(
        variable_name="delay_dry",
        display_name="Delay M-Beads dry",
        description="Delay protocol for drying of magnetic beads.",
        default=10,
        minimum=1,
        maximum=30,
        unit='min'
    )
    parameters.add_int(
        variable_name="delay_separate_final",
        display_name="Delay separate long",
        description="Delay protocol for separation of M-beads after adding elution buffer.",
        default=5,
        minimum=1,
        maximum=5,
        unit='min'
    )

def delay(protocol, delay_min, debug):
    if not debug:
        protocol.comment('Debugging mode: delay skipped.')
        protocol.delay(minutes=delay_min)
    protocol.comment(f'Delay protocol by {delay_min} min.')

class Task:
    def __init__(self, time_to_execute, action, sample_id):
        self.time_to_execute = time_to_execute  # Scheduled execution time (in seconds)
        self.action = action                    # Action to perform ('add_neutralization')
        self.sample_id = sample_id              # Identifier for the sample

    def __lt__(self, other):
        return self.time_to_execute < other.time_to_execute  # For heap ordering

def run(protocol: protocol_api.ProtocolContext):
    #variables not meant to be changed by user
    BOTTOM_DIST_DEFAULT_MM = 1
    BOTTOM_DIST_FAR_MM = 3
    
    
    # slots
    magnetic_module_slot = 1
    eppiracks_samples_slots = [2, 3]
    
    tipracks_300_slots = [5, 8, 11]
    tipracks_1000_slots = [4, 7, 10]
    
    tuberack_reagents_slot = 6
    
    # tipracks
    tips_300_per_sample = 8
    tips_1000_per_sample = 11
    
    tips_300_per_run = protocol.params.sample_count * tips_300_per_sample
    tips_1000_per_run = protocol.params.sample_count * tips_1000_per_sample
    
    tipracks_300_count = math.ceil(tips_300_per_run / 96)
    tipracks_1000_count = math.ceil(tips_1000_per_run / 96)
    
    tipracks_300 = [protocol.load_labware(load_name='opentrons_96_tiprack_300ul', location=slot) for slot in tipracks_300_slots[:tipracks_300_count]]
    tipracks_1000 = [protocol.load_labware(load_name='opentrons_96_tiprack_1000ul', location=slot) for slot in tipracks_1000_slots[:tipracks_1000_count]]
    
    # pipettes
    p300_single = protocol.load_instrument(instrument_name='p300_single', mount='right', tip_racks=tipracks_300)
    p1000_single = protocol.load_instrument(instrument_name='p1000_single', mount='left', tip_racks=tipracks_1000)
    
    # reagents
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
    
    # labware reagents
    tuberack_reagents = protocol.load_labware(load_name='opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical', location=tuberack_reagents_slot)
    
    tuberack_reagents[reagent_A1_slot].load_liquid(liquid=reagent_A1_liquid, volume=math.ceil((90 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_A2_slot].load_liquid(liquid=reagent_A2_liquid, volume=math.ceil((120 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_S3_slot].load_liquid(liquid=reagent_S3_liquid, volume=math.ceil((120 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_PAB_slot].load_liquid(liquid=reagent_PAB_liquid, volume=math.ceil((390 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_ERB_slot].load_liquid(liquid=reagent_ERB_liquid, volume=math.ceil((1800 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_AQ_slot].load_liquid(liquid=reagent_AQ_liquid, volume=math.ceil((1800 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_AE_slot].load_liquid(liquid=reagent_AE_liquid, volume=math.ceil((100 * protocol.params.sample_count * 1.1) / 10) * 10)
    tuberack_reagents[reagent_c_beads_slot].load_liquid(liquid=c_beads_liquid, volume=35 * protocol.params.sample_count * 1.1)
    tuberack_reagents[reagent_m_beads_slot].load_liquid(liquid=m_beads_liquid, volume=20 * protocol.params.sample_count * 1.1)
    
    # samples
    # unpurified_sample_eppi_rack_slot, unpurified_sample_slot,
    # purified_sample_eppi_rack_slot, purified_sample_slot,
    # processing_slot_1, processing_slot_2
    samples_infos = {
        0: [0, 'A1', 0, 'C1', 'A1', 'C1'],
        1: [0, 'A2', 0, 'C2', 'A2', 'C2'],
        2: [0, 'A3', 0, 'C3', 'A3', 'C3'],
        3: [0, 'A4', 0, 'C4', 'A4', 'C4'],
        4: [0, 'A5', 0, 'C5', 'A5', 'C5'],
        5: [0, 'A6', 0, 'C6', 'A6', 'C6'],
        6: [0, 'B1', 0, 'D1', 'A7', 'C7'],
        7: [0, 'B2', 0, 'D2', 'A8', 'C8'],
        8: [0, 'B3', 0, 'D3', 'A9', 'C9'],
        9: [0, 'B4', 0, 'D4', 'A10', 'C10'],
        10: [0, 'B5', 0, 'D5', 'A11', 'C11'],
        11: [0, 'B6', 0, 'D6', 'A12', 'C12'],
        12: [1, 'A1', 1, 'C1', 'B1', 'D1'],
        13: [1, 'A2', 1, 'C2', 'B2', 'D2'],
        14: [1, 'A3', 1, 'C3', 'B3', 'D3'],
        15: [1, 'A4', 1, 'C4', 'B4', 'D4'],
        16: [1, 'A5', 1, 'C5', 'B5', 'D5'],
        17: [1, 'A6', 1, 'C6', 'B6', 'D6'],
        18: [1, 'B1', 1, 'D1', 'B7', 'D7'],
        19: [1, 'B2', 1, 'D2', 'B8', 'D8'],
        20: [1, 'B3', 1, 'D3', 'B9', 'D9'],
        21: [1, 'B4', 1, 'D4', 'B10', 'D10'],
        22: [1, 'B5', 1, 'D5', 'B11', 'D11'],
        23: [1, 'B6', 1, 'D6', 'B12', 'D12']
    }

    # labware samples
    eppiracks_count = math.ceil(protocol.params.sample_count / 12)
    eppiracks_samples = [protocol.load_labware(load_name='opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', location=slot) for slot in eppiracks_samples_slots[:eppiracks_count]]
    
    for i in range(protocol.params.sample_count):
        bac_culture_sample_liquid = protocol.define_liquid(name="Bac. Pellet S" + str(i+1), description="Bacterial Culture Pellet S" + str(i+1), display_color="#00FF00")
        purified_plasmid_liquid = protocol.define_liquid(name="Purified Plasmid S" + str(i+1), description="Purified Plasmid S" + str(i+1), display_color="#0000FF")
        
        eppiracks_samples[samples_infos[i][0]][samples_infos[i][1]].load_liquid(liquid=bac_culture_sample_liquid, volume=20)
        eppiracks_samples[samples_infos[i][2]][samples_infos[i][3]].load_liquid(liquid=purified_plasmid_liquid, volume=0)

    # modules
    # magdeck = Magnetic Module GEN1
    mag_module = protocol.load_module('magdeck', magnetic_module_slot)
    mag_plate_96well = mag_module.load_labware('nest_96_wellplate_2ml_deep')
    
    #PROTOCOL START
    # Step 1: Add 90uL of reagent A1 to bacterial pellet, resuspend and transfer to magdeck 96well plate
    for i in range(protocol.params.sample_count):
        p300_single.transfer(90, tuberack_reagents[reagent_A1_slot].bottom(BOTTOM_DIST_DEFAULT_MM), eppiracks_samples[samples_infos[i][0]][samples_infos[i][1]].bottom(BOTTOM_DIST_DEFAULT_MM),
                             mix_after=(protocol.params.mix_times_resuspend_culture, 50), blow_out=True, blowout_location='destination_well')
        #transfer more volume than 90uL to account for extra volume of pellet
        p300_single.transfer(110, eppiracks_samples[samples_infos[i][0]][samples_infos[i][1]].bottom(BOTTOM_DIST_DEFAULT_MM), mag_plate_96well[samples_infos[i][4]].bottom(BOTTOM_DIST_DEFAULT_MM),
                             blow_out=True, blowout_location='destination_well')

    task_queue = []
    start_time = time.time()
    for i in range(protocol.params.sample_count):
        while task_queue and task_queue[0].time_to_execute < (time.time() - start_time):
            task = heapq.heappop(task_queue)
            if task.action == 'add_neutralization':
                # Step 3: Add 120uL neutralization buffer and mix
                p300_single.transfer(120, tuberack_reagents[reagent_S3_slot].bottom(BOTTOM_DIST_DEFAULT_MM), mag_plate_96well[samples_infos[task.sample_id][4]].bottom(BOTTOM_DIST_DEFAULT_MM),
                                     mix_after=(protocol.params.mix_times_default, 150), blow_out=True, blowout_location='destination_well')
        
        # Step 2: Add 120uL lysis buffer and mix, schedule neutralization step at delay_lysis time in the future (this is time critical)
        p300_single.transfer(120, tuberack_reagents[reagent_A2_slot].bottom(BOTTOM_DIST_DEFAULT_MM), mag_plate_96well[samples_infos[i][4]].bottom(BOTTOM_DIST_DEFAULT_MM),
                             mix_after=(protocol.params.mix_times_default, 100), blow_out=True, blowout_location='destination_well')
        neutralization_time = time.time() - start_time + protocol.params.delay_lysis*60 # Schedule after delay_lysis minutes 
        heapq.heappush(task_queue, Task(neutralization_time, 'add_neutralization', i))

    # process pending neutralization steps
    while task_queue:
        if task_queue[0].time_to_execute > (time.time() - start_time):
            delay(protocol, math.ceil(task_queue[0].time_to_execute - (time.time() - start_time))/60, protocol.params.debug)
            
        task = heapq.heappop(task_queue)
        if task.action == 'add_neutralization':
            # Step 3: Add 120uL neutralization buffer and mix
            p300_single.transfer(120, tuberack_reagents[reagent_S3_slot].bottom(BOTTOM_DIST_DEFAULT_MM), mag_plate_96well[samples_infos[task.sample_id][4]].bottom(BOTTOM_DIST_DEFAULT_MM),
                                 mix_after=(protocol.params.mix_times_default, 150), blow_out=True, blowout_location='destination_well')

    # Step 4: Add 35uL NucleoMag Clearing Beads and mix
    for i in range(protocol.params.sample_count):
        p300_single.transfer(35, tuberack_reagents[reagent_c_beads_slot].bottom(BOTTOM_DIST_DEFAULT_MM), mag_plate_96well[samples_infos[i][4]].bottom(BOTTOM_DIST_DEFAULT_MM),
                             mix_before=(protocol.params.mix_times_beads, 20), mix_after=(protocol.params.mix_times_thorough, 200), blow_out=True, blowout_location='destination_well')
    # not specified in protocol but feels right
    delay(protocol, protocol.params.delay_cbeads_incubate, protocol.params.debug)

    # Step 5-6: Magnetic separation and transfer supernatant to new well
    mag_module.engage(height_from_base=protocol.params.engage_height)
    delay(protocol, protocol.params.delay_separate_default, protocol.params.debug)
    for i in range(protocol.params.sample_count):
        p1000_single.transfer(365, mag_plate_96well[samples_infos[i][4]].bottom(BOTTOM_DIST_DEFAULT_MM), mag_plate_96well[samples_infos[i][5]].bottom(BOTTOM_DIST_DEFAULT_MM), blow_out=True, blowout_location='destination_well')
    mag_module.disengage()

    # Step 7: Add 20uL of NucleoMag M-Beads and 390uL of reagent PAB and mix
    for i in range(protocol.params.sample_count):
        p300_single.transfer(20, tuberack_reagents[reagent_m_beads_slot].bottom(BOTTOM_DIST_DEFAULT_MM), mag_plate_96well[samples_infos[i][5]].bottom(BOTTOM_DIST_DEFAULT_MM),
                             mix_before=(protocol.params.mix_times_beads, 10), blow_out=True, blowout_location='destination_well')
        p1000_single.transfer(390, tuberack_reagents[reagent_PAB_slot].bottom(BOTTOM_DIST_FAR_MM), mag_plate_96well[samples_infos[i][5]].bottom(BOTTOM_DIST_DEFAULT_MM),
                              mix_before=(protocol.params.mix_times_default, 390), mix_after=(protocol.params.mix_times_thorough, 400), blow_out=True, blowout_location='destination_well')
    delay(protocol, protocol.params.delay_resuspend_thorough, protocol.params.debug)
    
    # Step 8: Magnetic separation and remove supernatant
    mag_module.engage(height_from_base=protocol.params.engage_height)
    delay(protocol, protocol.params.delay_separate_default, protocol.params.debug)
    for i in range(protocol.params.sample_count):
        p1000_single.pick_up_tip()
        p1000_single.aspirate(775, mag_plate_96well[samples_infos[i][5]].bottom(BOTTOM_DIST_DEFAULT_MM))  # waste
        p1000_single.drop_tip()
    mag_module.disengage()

    # Steps 9-11: Wash with 900uL of ERB and AQ reagent, mix, remove supernatant and repeat
    # Wash 1 and 2 with ERB; Wash 3 and 4 with AQ 
    for j in range(4):
        if j == 0 or j == 1:
            for i in range(protocol.params.sample_count):
                p1000_single.transfer(900, tuberack_reagents[reagent_ERB_slot].bottom(BOTTOM_DIST_FAR_MM), mag_plate_96well[samples_infos[i][5]].bottom(BOTTOM_DIST_DEFAULT_MM),
                                      mix_after=(protocol.params.mix_times_thorough, 500), blow_out=True, blowout_location='destination_well')
        elif j == 2 or j == 3:
            for i in range(protocol.params.sample_count):
                p1000_single.transfer(900, tuberack_reagents[reagent_AQ_slot].bottom(BOTTOM_DIST_FAR_MM), mag_plate_96well[samples_infos[i][5]].bottom(BOTTOM_DIST_DEFAULT_MM),
                                      mix_after=(protocol.params.mix_times_thorough, 500), blow_out=True, blowout_location='destination_well')
            
        delay(protocol, protocol.params.delay_resuspend_wash, protocol.params.debug)
        mag_module.engage(height_from_base=protocol.params.engage_height)
        delay(protocol, protocol.params.delay_separate_default, protocol.params.debug)
        for i in range(protocol.params.sample_count):
            p1000_single.pick_up_tip()
            p1000_single.aspirate(900, mag_plate_96well[samples_infos[i][5]].bottom(BOTTOM_DIST_DEFAULT_MM))  # waste
            p1000_single.drop_tip()
        mag_module.disengage()


    # Step 16: Let beads dry for 15min
    delay(protocol, protocol.params.delay_dry, protocol.params.debug)
    
    # Step 17: Add 100uL of reagent AE, mix, and resuspend
    for i in range(protocol.params.sample_count):
        p300_single.transfer(100, tuberack_reagents[reagent_AE_slot].bottom(BOTTOM_DIST_DEFAULT_MM), mag_plate_96well[samples_infos[i][5]].bottom(BOTTOM_DIST_DEFAULT_MM),
                             mix_after=(protocol.params.mix_times_thorough, 50), blow_out=True, blowout_location='destination_well')

    # Step 18: Magnetic separation and transfer supernatant to final eppi for retrieval
    mag_module.engage(height_from_base=protocol.params.engage_height)
    delay(protocol, protocol.params.delay_separate_final, protocol.params.debug)
    for i in range(protocol.params.sample_count):
        p300_single.transfer(100, mag_plate_96well[samples_infos[i][5]].bottom(BOTTOM_DIST_DEFAULT_MM), eppiracks_samples[samples_infos[i][2]][samples_infos[i][3]].bottom(BOTTOM_DIST_DEFAULT_MM),
                             blow_out=True, blowout_location='destination_well')
    mag_module.disengage()