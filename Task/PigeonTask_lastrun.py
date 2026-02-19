#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2024.2.2),
    on Thu Feb 19 09:40:05 2026
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

# --- Import packages ---
from psychopy import locale_setup
from psychopy import prefs
from psychopy import plugins
plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, hardware
from psychopy.tools import environmenttools
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER, priority)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

from psychopy.hardware import keyboard

# --- Setup global variables (available in all functions) ---
# create a device manager to handle hardware (keyboards, mice, mirophones, speakers, etc.)
deviceManager = hardware.DeviceManager()
# ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
# store info about the experiment session
psychopyVersion = '2024.2.2'
expName = 'pigeon_task'  # from the Builder filename that created this script
# information about this experiment
expInfo = {
    'participant': '',
    'session': '001',
    'date|hid': data.getDateStr(),
    'expName|hid': expName,
    'psychopyVersion|hid': psychopyVersion,
}

# --- Define some variables which will change depending on pilot mode ---
'''
To run in pilot mode, either use the run/pilot toggle in Builder, Coder and Runner, 
or run the experiment with `--pilot` as an argument. To change what pilot 
#mode does, check out the 'Pilot mode' tab in preferences.
'''
# work out from system args whether we are running in pilot mode
PILOTING = core.setPilotModeFromArgs()
# start off with values from experiment settings
_fullScr = False
_winSize = [1440, 900]
# if in pilot mode, apply overrides according to preferences
if PILOTING:
    # force windowed mode
    if prefs.piloting['forceWindowed']:
        _fullScr = False
        # set window size
        _winSize = prefs.piloting['forcedWindowSize']

def showExpInfoDlg(expInfo):
    """
    Show participant info dialog.
    Parameters
    ==========
    expInfo : dict
        Information about this experiment.
    
    Returns
    ==========
    dict
        Information about this experiment.
    """
    # show participant info dialog
    dlg = gui.DlgFromDict(
        dictionary=expInfo, sortKeys=False, title=expName, alwaysOnTop=True
    )
    if dlg.OK == False:
        core.quit()  # user pressed cancel
    # return expInfo
    return expInfo


def setupData(expInfo, dataDir=None):
    """
    Make an ExperimentHandler to handle trials and saving.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    dataDir : Path, str or None
        Folder to save the data to, leave as None to create a folder in the current directory.    
    Returns
    ==========
    psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    # remove dialog-specific syntax from expInfo
    for key, val in expInfo.copy().items():
        newKey, _ = data.utils.parsePipeSyntax(key)
        expInfo[newKey] = expInfo.pop(key)
    
    # data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
    if dataDir is None:
        dataDir = _thisDir
    filename = u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])
    # make sure filename is relative to dataDir
    if os.path.isabs(filename):
        dataDir = os.path.commonprefix([dataDir, filename])
        filename = os.path.relpath(filename, dataDir)
    
    # an ExperimentHandler isn't essential but helps with data saving
    thisExp = data.ExperimentHandler(
        name=expName, version='',
        extraInfo=expInfo, runtimeInfo=None,
        originPath='/Users/jigold/GoldWorks/Mirror_jigold/Manuscripts/2023_Pigeon/Task/pigeontask_v7/PigeonTask_lastrun.py',
        savePickle=True, saveWideText=True,
        dataFileName=dataDir + os.sep + filename, sortColumns='time'
    )
    thisExp.setPriority('thisRow.t', priority.CRITICAL)
    thisExp.setPriority('expName', priority.LOW)
    # return experiment handler
    return thisExp


def setupLogging(filename):
    """
    Setup a log file and tell it what level to log at.
    
    Parameters
    ==========
    filename : str or pathlib.Path
        Filename to save log file and data files as, doesn't need an extension.
    
    Returns
    ==========
    psychopy.logging.LogFile
        Text stream to receive inputs from the logging system.
    """
    # set how much information should be printed to the console / app
    if PILOTING:
        logging.console.setLevel(
            prefs.piloting['pilotConsoleLoggingLevel']
        )
    else:
        logging.console.setLevel('warning')
    # save a log file for detail verbose info
    logFile = logging.LogFile(filename+'.log')
    if PILOTING:
        logFile.setLevel(
            prefs.piloting['pilotLoggingLevel']
        )
    else:
        logFile.setLevel(
            logging.getLevel('exp')
        )
    
    return logFile


def setupWindow(expInfo=None, win=None):
    """
    Setup the Window
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    win : psychopy.visual.Window
        Window to setup - leave as None to create a new window.
    
    Returns
    ==========
    psychopy.visual.Window
        Window in which to run this experiment.
    """
    if win is None:
        # if not given a window to setup, make one
        win = visual.Window(
            size=_winSize, fullscr=_fullScr, screen=0,
            winType='pyglet', allowStencil=True,
            monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
            backgroundImage='', backgroundFit='none',
            blendMode='avg', useFBO=True,
            units='height', 
            checkTiming=False  # we're going to do this ourselves in a moment
        )
    else:
        # if we have a window, just set the attributes which are safe to set
        win.color = [0,0,0]
        win.colorSpace = 'rgb'
        win.backgroundImage = ''
        win.backgroundFit = 'none'
        win.units = 'height'
    if expInfo is not None:
        # get/measure frame rate if not already in expInfo
        if win._monitorFrameRate is None:
            win._monitorFrameRate = win.getActualFrameRate(infoMsg='Attempting to measure frame rate of screen, please wait...')
        expInfo['frameRate'] = win._monitorFrameRate
    win.mouseVisible = True
    win.hideMessage()
    # show a visual indicator if we're in piloting mode
    if PILOTING and prefs.piloting['showPilotingIndicator']:
        win.showPilotingIndicator()
    
    return win


def setupDevices(expInfo, thisExp, win):
    """
    Setup whatever devices are available (mouse, keyboard, speaker, eyetracker, etc.) and add them to 
    the device manager (deviceManager)
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window in which to run this experiment.
    Returns
    ==========
    bool
        True if completed successfully.
    """
    # --- Setup input devices ---
    ioConfig = {}
    ioSession = ioServer = eyetracker = None
    
    # store ioServer object in the device manager
    deviceManager.ioServer = ioServer
    
    # create a default keyboard (e.g. to check for escape)
    if deviceManager.getDevice('defaultKeyboard') is None:
        deviceManager.addDevice(
            deviceClass='keyboard', deviceName='defaultKeyboard', backend='ptb'
        )
    if deviceManager.getDevice('sm_kbd_continue') is None:
        # initialise sm_kbd_continue
        sm_kbd_continue = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='sm_kbd_continue',
        )
    if deviceManager.getDevice('otr_kbd_response') is None:
        # initialise otr_kbd_response
        otr_kbd_response = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='otr_kbd_response',
        )
    if deviceManager.getDevice('pts_kbd_update') is None:
        # initialise pts_kbd_update
        pts_kbd_update = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='pts_kbd_update',
        )
    if deviceManager.getDevice('ptr_kbd_response') is None:
        # initialise ptr_kbd_response
        ptr_kbd_response = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='ptr_kbd_response',
        )
    if deviceManager.getDevice('sfm_kbd_continue') is None:
        # initialise sfm_kbd_continue
        sfm_kbd_continue = deviceManager.addDevice(
            deviceClass='keyboard',
            deviceName='sfm_kbd_continue',
        )
    # return True if completed successfully
    return True

def pauseExperiment(thisExp, win=None, timers=[], playbackComponents=[]):
    """
    Pause this experiment, preventing the flow from advancing to the next routine until resumed.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window for this experiment.
    timers : list, tuple
        List of timers to reset once pausing is finished.
    playbackComponents : list, tuple
        List of any components with a `pause` method which need to be paused.
    """
    # if we are not paused, do nothing
    if thisExp.status != PAUSED:
        return
    
    # start a timer to figure out how long we're paused for
    pauseTimer = core.Clock()
    # pause any playback components
    for comp in playbackComponents:
        comp.pause()
    # make sure we have a keyboard
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        defaultKeyboard = deviceManager.addKeyboard(
            deviceClass='keyboard',
            deviceName='defaultKeyboard',
            backend='PsychToolbox',
        )
    # run a while loop while we wait to unpause
    while thisExp.status == PAUSED:
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=['escape']):
            endExperiment(thisExp, win=win)
        # sleep 1ms so other threads can execute
        clock.time.sleep(0.001)
    # if stop was requested while paused, quit
    if thisExp.status == FINISHED:
        endExperiment(thisExp, win=win)
    # resume any playback components
    for comp in playbackComponents:
        comp.play()
    # reset any timers
    for timer in timers:
        timer.addTime(-pauseTimer.getTime())


def run(expInfo, thisExp, win, globalClock=None, thisSession=None):
    """
    Run the experiment flow.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    psychopy.visual.Window
        Window in which to run this experiment.
    globalClock : psychopy.core.clock.Clock or None
        Clock to get global time from - supply None to make a new one.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    # mark experiment as started
    thisExp.status = STARTED
    # make sure variables created by exec are available globally
    exec = environmenttools.setExecEnvironment(globals())
    # get device handles from dict of input devices
    ioServer = deviceManager.ioServer
    # get/create a default keyboard (e.g. to check for escape)
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        deviceManager.addDevice(
            deviceClass='keyboard', deviceName='defaultKeyboard', backend='PsychToolbox'
        )
    eyetracker = deviceManager.getDevice('eyetracker')
    # make sure we're running in the directory for this experiment
    os.chdir(_thisDir)
    # get filename from ExperimentHandler for convenience
    filename = thisExp.dataFileName
    frameTolerance = 0.001  # how close to onset before 'same' frame
    endExpNow = False  # flag for 'escape' or other condition => quit the exp
    # get frame duration from frame rate in expInfo
    if 'frameRate' in expInfo and expInfo['frameRate'] is not None:
        frameDur = 1.0 / round(expInfo['frameRate'])
    else:
        frameDur = 1.0 / 60.0  # could not measure, so guess
    
    # Start Code - component code to be run after the window creation
    
    # --- Initialize components for Routine "consent" ---
    consent_txt = visual.TextBox2(
         win, text='This is an academic research project conducted through the University of Pennsylvania. In this game, a pigeon will wander from the center of the screen to one of two seed piles on the left or right. Your goal is either to guess the direction the pigeon will ultimately choose or set a cut-off for the pigeon’s movement to make a choice based on their current location.\n\nThe estimated total time is about 30 minutes.\n\nYou must be at least 18 years old to participate. Your participation in this research is voluntary, which means you can choose whether or not to participate without negative consequences. Your anonymity is assured: the researchers who have requested your participation will not receive any personal information about you except your previously provided Prolific demographic data such as gender, ethnicity, and age. The de-identified data may be stored and distributed for future research studies without additional informed consent.\n\nIf you have questions about this research, please contact Alice Dallstream by emailing adalls@pennmedicine.upenn.edu or Josh Gold by emailing jigold@pennmedicine.upenn.edu.\n\nIf you have questions, concerns, or complaints regarding your participation in this research study, or if you have any questions about your rights as a research subject and you cannot reach a member of the study team, you may contact the Office of Regulatory Affairs at the University of Pennsylvania by calling (215) 898-2614 or emailing irb@pobox.upenn.edu.\n', placeholder='Type here...', font='Open Sans',
         ori=0.0, pos=(0, 0), draggable=False,      letterHeight=0.023,
         size=(0.75, 0.68), borderWidth=2.0,
         color='white', colorSpace='rgb',
         opacity=None,
         bold=False, italic=False,
         lineSpacing=1.0, speechPoint=None,
         padding=0.0, alignment='center',
         anchor='center', overflow='visible',
         fillColor=[-0.1765, -0.1765, -0.1765], borderColor=None,
         flipHoriz=False, flipVert=False, languageStyle='LTR',
         editable=False,
         name='consent_txt',
         depth=0, autoLog=True,
    )
    consent_box = visual.Rect(
        win=win, name='consent_box',
        width=(0.05, 0.05)[0], height=(0.05, 0.05)[1],
        ori=0.0, pos=(-0.4, -0.4), draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor='white',
        opacity=None, depth=-1.0, interpolate=True)
    consent_mouse = event.Mouse(win=win)
    x, y = [None, None]
    consent_mouse.mouseClock = core.Clock()
    consent_acknowledge = visual.TextBox2(
         win, text='I have read the information above and agree to take part in this study. I understand that I may withdraw my consent at any time before I complete the tasks.', placeholder='Type here...', font='Open Sans',
         ori=0.0, pos=(-0.36, -0.42), draggable=False,      letterHeight=0.023,
         size=(0.71, 0.08), borderWidth=2.0,
         color='white', colorSpace='rgb',
         opacity=None,
         bold=False, italic=False,
         lineSpacing=1.0, speechPoint=None,
         padding=0.0, alignment='center',
         anchor='center-left', overflow='visible',
         fillColor=[-0.1765, -0.1765, -0.1765], borderColor=None,
         flipHoriz=False, flipVert=False, languageStyle='LTR',
         editable=False,
         name='consent_acknowledge',
         depth=-3, autoLog=True,
    )
    
    # --- Initialize components for Routine "experiment_setup" ---
    # Run 'Begin Experiment' code from es_code_setup
    # Global CONSTANTS and variables
    
    # Local variables needed in JS:
    still_going = True
    choice = ''
    
    # Colors
    # For JS:
    #white = new util.Color([1, 1, 1]);
    #grey = new util.Color([0.2, 0.2, 0.2]);
    #yellow = new util.Color([1, 1, 0]);
    #green = new util.Color([-1, 0, -1]);
    #lightgreen = new util.Color([0, 1, 0]);
    #black = new util.Color([-1, -1, -1]);
    #red = new util.Color([1, 0, 0]);
    white=[1,1,1]
    grey=[0.2,0.2,0.2]
    yellow=[1,1,0]
    green=[-1,0,-1]
    lightgreen=[0,1,0]
    black=[-1,-1,-1]
    red=[1,0,0]
    
    # Screen
    ORIGIN_X = 0
    ORIGIN_Y = 0
    SEEDS_WIDTH = 0.1
    SEEDS_HEIGHT = 0.1
    EDGE_DISTANCE = 0.75
    
    # Messages
    INSTRUCTION_X = ORIGIN_X
    INSTRUCTION_Y = 0.45
    INSTRUCTION_COLOR = yellow
    CHOICE_STRING = 'Press Left or Right arrow to choose'
    CONTINUE_STRING = 'Press Spacebar or any arrow to continue'
    PREDEFINED_CHOICE_STRING='Press Left/Up for small/large increases, Right/Down for small/large decreases,\nSpacebar to start the trial.'
    PREDEFINED_CHOICE_COMMIT_STRING='Press Left/Up for small/large increases, Right/Down for small/large decreases,\nc to commit or Spacebar to start the trial.'
    PREDEFINED_NO_CHOICE_STRING='Petey failed to reach the bound\nConsider moving it closer.'
    PREDEFINED_ABORT_INSTRUCTION='Press Spacebar or any arrow to abort trial'
    PREDEFINED_ABORT_STRING='Aborted early'
    CORRECT_STRING='Correct'
    ERROR_STRING='Error'
    MESSAGE_X = ORIGIN_X
    MESSAGE_Y = 0.3
    instruction_show = True
    instruction_string = CONTINUE_STRING
    message_show = True
    message_string = ''
    steps_string = ' '
    coins_string = ' '
    
    # Status symbols
    LINE_HEIGHT = 0.1
    STATUS_BAR_HEIGHT = 0.2
    STATUS_BAR_WIDTH = 0.1
    STATUS_BAR_BASE_Y = ORIGIN_Y-0.4
    STATUS_BAR_APEX_Y = STATUS_BAR_BASE_Y+STATUS_BAR_HEIGHT
    STATUS_BAR_BRACKET_WIDTH =  0.35
    ABSCISSA_Y = ORIGIN_Y-0.05
    ABSCISSA_WIDTH = EDGE_DISTANCE*2
    ORDINATE_Y = ORIGIN_Y+0.05
    ORDINATE_WIDTH = 0.2
    ORDINATE_HEIGHT = 0.01
    COINS_BAR_X = -0.1
    COINS_TEXT_X = -0.1
    COINS_TEXT_Y = -0.37
    COINS_LABEL_TEXT_X = -0.28
    COINS_LABEL_TEXT_Y = -0.2
    COINS_PROPERTIES_TEXT_X = -0.55
    COINS_PROPERTIES_TEXT_Y = -0.3
    STEPS_BAR_X = 0.1
    STEPS_TEXT_X = 0.1
    STEPS_TEXT_Y = -0.37
    STEPS_LABEL_TEXT_X = 0.3
    STEPS_LABEL_TEXT_Y = -0.2
    STEPS_PROPERTIES_TEXT_X = 0.55
    STEPS_PROPERTIES_TEXT_Y = -0.3
    steps_bar_height = 0
    steps_bar_center = 0
    coins_bar_height = 0
    coins_bar_center = 0
    
    # Predefined bound
    PREDEFINED_BOUND_COLOR = green
    PREDEFINED_BOUND_WIDTH = 0.2
    PREDEFINED_BOUND_HEIGHT = 0.01
    PREDEFINED_BOUND_Y = ORIGIN_Y+0.05
    PREDEFINED_BOUND_SMALL_DELTA = 0.003
    PREDEFINED_BOUND_LARGE_DELTA = 0.05
    PREDEFINED_BOUND_X_DEFAULT = 0.05
    PREDEFINED_BOUND_MARKER_SZ = 0.01
    PREDEFINED_BOUND_MARKER_Y = 0.16
    PREDEFINED_BOUND_MARKER_COLOR = lightgreen
    predefined_bound_x = PREDEFINED_BOUND_X_DEFAULT
    predefined_bound_x_previous = predefined_bound_x
    committed_to_bound = False
    predefined_bound_min_bonus = 0
    predefined_bound_max_bonus = 0
    min_steps_to_commit = 100
    
    # Stimulus display
    seeds_show = False
    pigeon_show = False
    status_show = False
    box_show = False
    
    # Timing/counting
    UPDATE_INTERVAL = 0.2
    PREDEFINED_UPDATE_INTERVAL = 0.1
    block_number = 0
    trial_number = 0
    trial_start = True
    bonus_count = 0
    
    # Task conditions
    task_type = 'online'
    STEPS_MAX_PER_TRIAL = 50
    coins_gained_per_correct = 1
    coins_lost_per_error = 1
    coins_paid_to_start_trial = 0
    coins_count = 2
    coins_max = 8
    steps_lost_per_error = 0
    steps_taken_to_start_trial = 1
    steps_count = 0
    steps_max = 200
    step_mean_val = 0.02
    step_std_val = 0.1
    last_coins_count = 0
    last_steps_count = 0
    
    # Pigeon position
    pigeon_flip = 0.1
    pigeon_true_x = 0
    pigeon_shown_x = pigeon_true_x
    pigeon_steps = []
    
    # Function to move pigeon, depends on random variables
    #   Note -- if you change this code, set code type to 
    #   "Auto->JS" first to make sure the JS code gets updated,
    #   then set code type back to "both" and add the following
    #   line at the top of this function:
    #    global pigeon_flip, pigeon_true_x, pigeon_shown_x
    def move_pigeon_function(step_mean, step_std):
        global pigeon_flip, pigeon_true_x, pigeon_shown_x
        
        for i in range(10):    
            u = 2*random()-1
            v = 2*random()-1
            r = u*u + v*v
            if r!=0 and r<1:
                break
        new_step = step_mean + step_std*u*sqrt(-2*log(r)/r)
        if new_step > 0:
            pigeon_flip = -0.1
        else:
            pigeon_flip = 0.1
    
        pigeon_true_x = pigeon_true_x + new_step
        if pigeon_true_x >= EDGE_DISTANCE:
            pigeon_shown_x = EDGE_DISTANCE
        elif pigeon_true_x <= -EDGE_DISTANCE:
            pigeon_shown_x = -EDGE_DISTANCE
        else:
            pigeon_shown_x = pigeon_true_x
    
    # Function to update status text/graphics
    #   Note -- if you change this code, set code type to 
    #   "Auto->JS" first to make sure the JS code gets updated,
    #   then set code type back to "both" and add the following
    #   line at the top of this function:
    #    global coins_string, steps_string, coins_bar_height, coins_bar_center, steps_bar_height, steps_bar_center
    def update_status_function(update_text, update_coins_bar, update_steps_bar):
        global coins_string, steps_string, coins_bar_height, coins_bar_center, steps_bar_height, steps_bar_center
    
        # Conditionally update strings
        if update_text:
            coins_string = 'Coins gained per correct choice = ' + str(coins_gained_per_correct) + '\nCoins lost per error choice = ' + str(coins_lost_per_error) + '\nCoins paid to start each trial = ' + str(coins_paid_to_start_trial) + '\nTarget number of coins for bonus = ' + str(coins_max)
            steps_string = 'Total steps = ' + str(steps_max) + '\nSteps lost per error choice = ' + str(steps_lost_per_error) + '\nSteps taken to start each trial = ' + str(steps_taken_to_start_trial)
    
        # Conditionally update coins bar
        if update_coins_bar:
            if coins_count <= coins_max:
                coins_bar_height = coins_count/coins_max*STATUS_BAR_HEIGHT
            else:
                coins_bar_height = STATUS_BAR_HEIGHT+0.02
            coins_bar_center = STATUS_BAR_BASE_Y+coins_bar_height/2
    
        # Conditionally update steps bar
        if update_steps_bar:
            steps_bar_height = (steps_max-steps_count)/steps_max*STATUS_BAR_HEIGHT
            steps_bar_center = STATUS_BAR_BASE_Y+steps_bar_height/2
    
    # Function to show feedback after each trial
    #   Note -- if you change this code, set code type to 
    #   "Auto->JS" first to make sure the JS code gets updated,
    #   then set code type back to "both" and add the following
    #   line at the top of this function:
    #   global last_coins_count, last_steps_count, message_string
    def set_trial_end_string_function(feedback_string):
        global last_coins_count, last_steps_count, message_string
    
        coins_change = coins_count - last_coins_count
        steps_change = steps_count - last_steps_count
        if coins_change > 0:
            coins_feedback_string = str(coins_change) + ' coins gained'
        elif coins_change < 0:
            coins_feedback_string = str(coins_change) + ' coins lost'
        else:
            coins_feedback_string = 'no coins gained or lost'
        
        message_string = feedback_string + '\n\n' + coins_feedback_string + '\n\n' + str(steps_change) + ' steps used this trial'
        last_coins_count = coins_count
        last_steps_count = steps_count
        
    # Function to generate standard string at the end of each block
    #   Note -- if you change this code, set code type to 
    #   "Auto->JS" first to make sure the JS code gets updated,
    #   then set code type back to "both" and add the following
    #   line at the top of this function:
    #    global message_string, bonus_count, coins_count, coins_bar_height, coins_bar_center
    def set_block_end_string_function(): 
        global message_string, bonus_count, coins_count, coins_bar_height, coins_bar_center
    
        if (task_type == 'predefined_commit' and predefined_bound_x >= predefined_bound_min_bonus and predefined_bound_x <= predefined_bound_max_bonus) or (task_type != 'predefined_commit' and coins_count >= coins_max):                
            message_string='Block finished.\n\nYou earned a bonus. Great job!'
            bonus_count += 1
            coins_bar_height = STATUS_BAR_HEIGHT
            coins_bar_center = STATUS_BAR_BASE_Y+coins_bar_height/2
            coins_count = 'Bonus!'
        else:
            message_string='Block finished.\n\nGreat job!'
            coins_bar_height = 0
            coins_bar_center = STATUS_BAR_BASE_Y
            coins_count = 'No bonus'
    
    move_pigeon = move_pigeon_function
    update_status = update_status_function
    set_trial_end_string = set_trial_end_string_function
    set_block_end_string = set_block_end_string_function
    update_status(True, True, True)
    
    # --- Initialize components for Routine "general_instructions_setup" ---
    
    # --- Initialize components for Routine "show_message" ---
    sm_kbd_continue = keyboard.Keyboard(deviceName='sm_kbd_continue')
    # Run 'Begin Experiment' code from sm_code_update
    box_width = 0
    box_height = 0
    box_x = 0
    box_y = 0
    sm_polygon_box = visual.Rect(
        win=win, name='sm_polygon_box',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-2.0, interpolate=True)
    sm_polygon_abscissa = visual.Line(
        win=win, name='sm_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    sm_polygon_ordinate = visual.Line(
        win=win, name='sm_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-4.0, interpolate=True)
    sm_polygon_petey_midpoint = visual.Line(
        win=win, name='sm_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-5.0, interpolate=True)
    sm_image_right_seeds = visual.ImageStim(
        win=win,
        name='sm_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-6.0)
    sm_image_left_seeds = visual.ImageStim(
        win=win,
        name='sm_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    sm_image_petey = visual.ImageStim(
        win=win,
        name='sm_image_petey', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-8.0)
    sm_polygon_coin_bar = visual.Rect(
        win=win, name='sm_polygon_coin_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-9.0, interpolate=True)
    sm_polygon_steps_bar = visual.Rect(
        win=win, name='sm_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-10.0, interpolate=True)
    sm_polygon_status_bar_apex = visual.Line(
        win=win, name='sm_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-11.0, interpolate=True)
    sm_polygon_status_bar_base = visual.Line(
        win=win, name='sm_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-12.0, interpolate=True)
    sm_text_coins = visual.TextStim(win=win, name='sm_text_coins',
        text='',
        font='Open Sans',
        pos=(COINS_TEXT_X,COINS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-13.0);
    sm_text_coins_label = visual.TextStim(win=win, name='sm_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-14.0);
    sm_text_coins_properties = visual.TextStim(win=win, name='sm_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-15.0);
    sm_text_steps = visual.TextStim(win=win, name='sm_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-16.0);
    sm_text_steps_label = visual.TextStim(win=win, name='sm_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-17.0);
    sm_text_steps_properties = visual.TextStim(win=win, name='sm_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-18.0);
    sm_text_message = visual.TextStim(win=win, name='sm_text_message',
        text='',
        font='Open Sans',
        pos=(MESSAGE_X,MESSAGE_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-19.0);
    sm_text_instruction = visual.TextStim(win=win, name='sm_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-20.0);
    
    # --- Initialize components for Routine "simple_demo_setup" ---
    
    # --- Initialize components for Routine "show_message" ---
    sm_kbd_continue = keyboard.Keyboard(deviceName='sm_kbd_continue')
    # Run 'Begin Experiment' code from sm_code_update
    box_width = 0
    box_height = 0
    box_x = 0
    box_y = 0
    sm_polygon_box = visual.Rect(
        win=win, name='sm_polygon_box',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-2.0, interpolate=True)
    sm_polygon_abscissa = visual.Line(
        win=win, name='sm_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    sm_polygon_ordinate = visual.Line(
        win=win, name='sm_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-4.0, interpolate=True)
    sm_polygon_petey_midpoint = visual.Line(
        win=win, name='sm_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-5.0, interpolate=True)
    sm_image_right_seeds = visual.ImageStim(
        win=win,
        name='sm_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-6.0)
    sm_image_left_seeds = visual.ImageStim(
        win=win,
        name='sm_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    sm_image_petey = visual.ImageStim(
        win=win,
        name='sm_image_petey', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-8.0)
    sm_polygon_coin_bar = visual.Rect(
        win=win, name='sm_polygon_coin_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-9.0, interpolate=True)
    sm_polygon_steps_bar = visual.Rect(
        win=win, name='sm_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-10.0, interpolate=True)
    sm_polygon_status_bar_apex = visual.Line(
        win=win, name='sm_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-11.0, interpolate=True)
    sm_polygon_status_bar_base = visual.Line(
        win=win, name='sm_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-12.0, interpolate=True)
    sm_text_coins = visual.TextStim(win=win, name='sm_text_coins',
        text='',
        font='Open Sans',
        pos=(COINS_TEXT_X,COINS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-13.0);
    sm_text_coins_label = visual.TextStim(win=win, name='sm_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-14.0);
    sm_text_coins_properties = visual.TextStim(win=win, name='sm_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-15.0);
    sm_text_steps = visual.TextStim(win=win, name='sm_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-16.0);
    sm_text_steps_label = visual.TextStim(win=win, name='sm_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-17.0);
    sm_text_steps_properties = visual.TextStim(win=win, name='sm_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-18.0);
    sm_text_message = visual.TextStim(win=win, name='sm_text_message',
        text='',
        font='Open Sans',
        pos=(MESSAGE_X,MESSAGE_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-19.0);
    sm_text_instruction = visual.TextStim(win=win, name='sm_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-20.0);
    
    # --- Initialize components for Routine "simple_demo" ---
    sd_polygon_abscissa = visual.Line(
        win=win, name='sd_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-1.0, interpolate=True)
    sd_image_left_seeds = visual.ImageStim(
        win=win,
        name='sd_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-2.0)
    sd_image_right_seeds = visual.ImageStim(
        win=win,
        name='sd_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-3.0)
    sd_image_pigeon = visual.ImageStim(
        win=win,
        name='sd_image_pigeon', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-4.0)
    
    # --- Initialize components for Routine "status_instruction_setup" ---
    
    # --- Initialize components for Routine "show_message" ---
    sm_kbd_continue = keyboard.Keyboard(deviceName='sm_kbd_continue')
    # Run 'Begin Experiment' code from sm_code_update
    box_width = 0
    box_height = 0
    box_x = 0
    box_y = 0
    sm_polygon_box = visual.Rect(
        win=win, name='sm_polygon_box',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-2.0, interpolate=True)
    sm_polygon_abscissa = visual.Line(
        win=win, name='sm_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    sm_polygon_ordinate = visual.Line(
        win=win, name='sm_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-4.0, interpolate=True)
    sm_polygon_petey_midpoint = visual.Line(
        win=win, name='sm_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-5.0, interpolate=True)
    sm_image_right_seeds = visual.ImageStim(
        win=win,
        name='sm_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-6.0)
    sm_image_left_seeds = visual.ImageStim(
        win=win,
        name='sm_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    sm_image_petey = visual.ImageStim(
        win=win,
        name='sm_image_petey', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-8.0)
    sm_polygon_coin_bar = visual.Rect(
        win=win, name='sm_polygon_coin_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-9.0, interpolate=True)
    sm_polygon_steps_bar = visual.Rect(
        win=win, name='sm_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-10.0, interpolate=True)
    sm_polygon_status_bar_apex = visual.Line(
        win=win, name='sm_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-11.0, interpolate=True)
    sm_polygon_status_bar_base = visual.Line(
        win=win, name='sm_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-12.0, interpolate=True)
    sm_text_coins = visual.TextStim(win=win, name='sm_text_coins',
        text='',
        font='Open Sans',
        pos=(COINS_TEXT_X,COINS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-13.0);
    sm_text_coins_label = visual.TextStim(win=win, name='sm_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-14.0);
    sm_text_coins_properties = visual.TextStim(win=win, name='sm_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-15.0);
    sm_text_steps = visual.TextStim(win=win, name='sm_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-16.0);
    sm_text_steps_label = visual.TextStim(win=win, name='sm_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-17.0);
    sm_text_steps_properties = visual.TextStim(win=win, name='sm_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-18.0);
    sm_text_message = visual.TextStim(win=win, name='sm_text_message',
        text='',
        font='Open Sans',
        pos=(MESSAGE_X,MESSAGE_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-19.0);
    sm_text_instruction = visual.TextStim(win=win, name='sm_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-20.0);
    
    # --- Initialize components for Routine "start_task" ---
    
    # --- Initialize components for Routine "show_message" ---
    sm_kbd_continue = keyboard.Keyboard(deviceName='sm_kbd_continue')
    # Run 'Begin Experiment' code from sm_code_update
    box_width = 0
    box_height = 0
    box_x = 0
    box_y = 0
    sm_polygon_box = visual.Rect(
        win=win, name='sm_polygon_box',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-2.0, interpolate=True)
    sm_polygon_abscissa = visual.Line(
        win=win, name='sm_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    sm_polygon_ordinate = visual.Line(
        win=win, name='sm_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-4.0, interpolate=True)
    sm_polygon_petey_midpoint = visual.Line(
        win=win, name='sm_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-5.0, interpolate=True)
    sm_image_right_seeds = visual.ImageStim(
        win=win,
        name='sm_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-6.0)
    sm_image_left_seeds = visual.ImageStim(
        win=win,
        name='sm_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    sm_image_petey = visual.ImageStim(
        win=win,
        name='sm_image_petey', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-8.0)
    sm_polygon_coin_bar = visual.Rect(
        win=win, name='sm_polygon_coin_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-9.0, interpolate=True)
    sm_polygon_steps_bar = visual.Rect(
        win=win, name='sm_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-10.0, interpolate=True)
    sm_polygon_status_bar_apex = visual.Line(
        win=win, name='sm_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-11.0, interpolate=True)
    sm_polygon_status_bar_base = visual.Line(
        win=win, name='sm_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-12.0, interpolate=True)
    sm_text_coins = visual.TextStim(win=win, name='sm_text_coins',
        text='',
        font='Open Sans',
        pos=(COINS_TEXT_X,COINS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-13.0);
    sm_text_coins_label = visual.TextStim(win=win, name='sm_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-14.0);
    sm_text_coins_properties = visual.TextStim(win=win, name='sm_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-15.0);
    sm_text_steps = visual.TextStim(win=win, name='sm_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-16.0);
    sm_text_steps_label = visual.TextStim(win=win, name='sm_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-17.0);
    sm_text_steps_properties = visual.TextStim(win=win, name='sm_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-18.0);
    sm_text_message = visual.TextStim(win=win, name='sm_text_message',
        text='',
        font='Open Sans',
        pos=(MESSAGE_X,MESSAGE_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-19.0);
    sm_text_instruction = visual.TextStim(win=win, name='sm_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-20.0);
    
    # --- Initialize components for Routine "setup_instructions_message" ---
    
    # --- Initialize components for Routine "show_message" ---
    sm_kbd_continue = keyboard.Keyboard(deviceName='sm_kbd_continue')
    # Run 'Begin Experiment' code from sm_code_update
    box_width = 0
    box_height = 0
    box_x = 0
    box_y = 0
    sm_polygon_box = visual.Rect(
        win=win, name='sm_polygon_box',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-2.0, interpolate=True)
    sm_polygon_abscissa = visual.Line(
        win=win, name='sm_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    sm_polygon_ordinate = visual.Line(
        win=win, name='sm_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-4.0, interpolate=True)
    sm_polygon_petey_midpoint = visual.Line(
        win=win, name='sm_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-5.0, interpolate=True)
    sm_image_right_seeds = visual.ImageStim(
        win=win,
        name='sm_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-6.0)
    sm_image_left_seeds = visual.ImageStim(
        win=win,
        name='sm_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    sm_image_petey = visual.ImageStim(
        win=win,
        name='sm_image_petey', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-8.0)
    sm_polygon_coin_bar = visual.Rect(
        win=win, name='sm_polygon_coin_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-9.0, interpolate=True)
    sm_polygon_steps_bar = visual.Rect(
        win=win, name='sm_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-10.0, interpolate=True)
    sm_polygon_status_bar_apex = visual.Line(
        win=win, name='sm_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-11.0, interpolate=True)
    sm_polygon_status_bar_base = visual.Line(
        win=win, name='sm_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-12.0, interpolate=True)
    sm_text_coins = visual.TextStim(win=win, name='sm_text_coins',
        text='',
        font='Open Sans',
        pos=(COINS_TEXT_X,COINS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-13.0);
    sm_text_coins_label = visual.TextStim(win=win, name='sm_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-14.0);
    sm_text_coins_properties = visual.TextStim(win=win, name='sm_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-15.0);
    sm_text_steps = visual.TextStim(win=win, name='sm_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-16.0);
    sm_text_steps_label = visual.TextStim(win=win, name='sm_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-17.0);
    sm_text_steps_properties = visual.TextStim(win=win, name='sm_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-18.0);
    sm_text_message = visual.TextStim(win=win, name='sm_text_message',
        text='',
        font='Open Sans',
        pos=(MESSAGE_X,MESSAGE_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-19.0);
    sm_text_instruction = visual.TextStim(win=win, name='sm_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-20.0);
    
    # --- Initialize components for Routine "online_task_run" ---
    otr_kbd_response = keyboard.Keyboard(deviceName='otr_kbd_response')
    otr_polygon_abscissa = visual.Line(
        win=win, name='otr_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-2.0, interpolate=True)
    otr_polygon_ordinate = visual.Line(
        win=win, name='otr_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    otr_polygon_petey_midpoint = visual.Line(
        win=win, name='otr_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor='red', fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-4.0, interpolate=True)
    otr_image_pigeon = visual.ImageStim(
        win=win,
        name='otr_image_pigeon', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-5.0)
    otr_image_right_seeds = visual.ImageStim(
        win=win,
        name='otr_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-6.0)
    otr_image_left_seeds = visual.ImageStim(
        win=win,
        name='otr_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    otr_polygon_coins_bar = visual.Rect(
        win=win, name='otr_polygon_coins_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor='white', fillColor='white',
        opacity=None, depth=-8.0, interpolate=True)
    otr_polygon_steps_bar = visual.Rect(
        win=win, name='otr_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor='white', fillColor='white',
        opacity=None, depth=-9.0, interpolate=True)
    otr_polygon_status_bar_apex = visual.Line(
        win=win, name='otr_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor='black', fillColor='black',
        opacity=None, depth=-10.0, interpolate=True)
    otr_polygon_status_bar_base = visual.Line(
        win=win, name='otr_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-11.0, interpolate=True)
    otr_text_coins = visual.TextStim(win=win, name='otr_text_coins',
        text='',
        font='Open Sans',
        pos=(COINS_TEXT_X,COINS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-12.0);
    otr_text_coins_label = visual.TextStim(win=win, name='otr_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-13.0);
    otr_text_coins_properties = visual.TextStim(win=win, name='otr_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-14.0);
    otr_text_steps = visual.TextStim(win=win, name='otr_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-15.0);
    otr_text_steps_label = visual.TextStim(win=win, name='otr_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X,STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-16.0);
    otr_text_steps_properties = visual.TextStim(win=win, name='otr_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-17.0);
    otr_text_instruction = visual.TextStim(win=win, name='otr_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-18.0);
    
    # --- Initialize components for Routine "show_message" ---
    sm_kbd_continue = keyboard.Keyboard(deviceName='sm_kbd_continue')
    # Run 'Begin Experiment' code from sm_code_update
    box_width = 0
    box_height = 0
    box_x = 0
    box_y = 0
    sm_polygon_box = visual.Rect(
        win=win, name='sm_polygon_box',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-2.0, interpolate=True)
    sm_polygon_abscissa = visual.Line(
        win=win, name='sm_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    sm_polygon_ordinate = visual.Line(
        win=win, name='sm_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-4.0, interpolate=True)
    sm_polygon_petey_midpoint = visual.Line(
        win=win, name='sm_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-5.0, interpolate=True)
    sm_image_right_seeds = visual.ImageStim(
        win=win,
        name='sm_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-6.0)
    sm_image_left_seeds = visual.ImageStim(
        win=win,
        name='sm_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    sm_image_petey = visual.ImageStim(
        win=win,
        name='sm_image_petey', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-8.0)
    sm_polygon_coin_bar = visual.Rect(
        win=win, name='sm_polygon_coin_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-9.0, interpolate=True)
    sm_polygon_steps_bar = visual.Rect(
        win=win, name='sm_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-10.0, interpolate=True)
    sm_polygon_status_bar_apex = visual.Line(
        win=win, name='sm_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-11.0, interpolate=True)
    sm_polygon_status_bar_base = visual.Line(
        win=win, name='sm_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-12.0, interpolate=True)
    sm_text_coins = visual.TextStim(win=win, name='sm_text_coins',
        text='',
        font='Open Sans',
        pos=(COINS_TEXT_X,COINS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-13.0);
    sm_text_coins_label = visual.TextStim(win=win, name='sm_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-14.0);
    sm_text_coins_properties = visual.TextStim(win=win, name='sm_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-15.0);
    sm_text_steps = visual.TextStim(win=win, name='sm_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-16.0);
    sm_text_steps_label = visual.TextStim(win=win, name='sm_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-17.0);
    sm_text_steps_properties = visual.TextStim(win=win, name='sm_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-18.0);
    sm_text_message = visual.TextStim(win=win, name='sm_text_message',
        text='',
        font='Open Sans',
        pos=(MESSAGE_X,MESSAGE_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-19.0);
    sm_text_instruction = visual.TextStim(win=win, name='sm_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-20.0);
    
    # --- Initialize components for Routine "predefined_task_setup" ---
    pts_kbd_update = keyboard.Keyboard(deviceName='pts_kbd_update')
    pts_polygon_abscissa = visual.Line(
        win=win, name='pts_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-2.0, interpolate=True)
    pts_polygon_ordinate = visual.Line(
        win=win, name='pts_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    pts_polygon_petey_midpoint = visual.Line(
        win=win, name='pts_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor='red', fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-4.0, interpolate=True)
    pts_polygon_left_bound = visual.Line(
        win=win, name='pts_polygon_left_bound',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor='white', fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-5.0, interpolate=True)
    pts_polygon_right_bound = visual.Line(
        win=win, name='pts_polygon_right_bound',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor=PREDEFINED_BOUND_COLOR, fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-6.0, interpolate=True)
    pts_polygon_stay_left = visual.Line(
        win=win, name='pts_polygon_stay_left',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=PREDEFINED_BOUND_MARKER_COLOR, fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-7.0, interpolate=True)
    pts_polygon_stay_right = visual.Line(
        win=win, name='pts_polygon_stay_right',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=PREDEFINED_BOUND_MARKER_COLOR, fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-8.0, interpolate=True)
    pts_image_pigeon = visual.ImageStim(
        win=win,
        name='pts_image_pigeon', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-9.0)
    pts_image_right_seeds = visual.ImageStim(
        win=win,
        name='pts_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-10.0)
    pts_image_left_seeds = visual.ImageStim(
        win=win,
        name='pts_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-11.0)
    pts_polygon_coins_bar = visual.Rect(
        win=win, name='pts_polygon_coins_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor='white', fillColor='white',
        opacity=None, depth=-12.0, interpolate=True)
    pts_polygon_steps_bar = visual.Rect(
        win=win, name='pts_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor='white', fillColor='white',
        opacity=None, depth=-13.0, interpolate=True)
    pts_polygon_status_bar_apex = visual.Line(
        win=win, name='pts_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor='black', fillColor='black',
        opacity=None, depth=-14.0, interpolate=True)
    pts_polygon_status_bar_base = visual.Line(
        win=win, name='pts_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-15.0, interpolate=True)
    pts_text_coins = visual.TextStim(win=win, name='pts_text_coins',
        text='',
        font='Open Sans',
        pos=[0,0], draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-16.0);
    pts_text_coins_label = visual.TextStim(win=win, name='pts_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-17.0);
    pts_text_coins_properties = visual.TextStim(win=win, name='pts_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-18.0);
    pts_text_steps = visual.TextStim(win=win, name='pts_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-19.0);
    pts_text_steps_label = visual.TextStim(win=win, name='pts_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X,STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-20.0);
    pts_text_steps_properties = visual.TextStim(win=win, name='pts_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-21.0);
    pts_text_instruction = visual.TextStim(win=win, name='pts_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-22.0);
    pts_text_hint = visual.TextStim(win=win, name='pts_text_hint',
        text='-> move closer to save steps <-\n<- move farther to increase accuracy ->',
        font='Open Sans',
        pos=[0,0], draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-23.0);
    pts_text_hint_2 = visual.TextStim(win=win, name='pts_text_hint_2',
        text='The light green tics indicate your last choice.',
        font='Open Sans',
        pos=[0,0], draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=[0.1294, 0.8667, 0.1294], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-24.0);
    
    # --- Initialize components for Routine "predefined_task_run" ---
    ptr_kbd_response = keyboard.Keyboard(deviceName='ptr_kbd_response')
    ptr_polygon_abscissa = visual.Line(
        win=win, name='ptr_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-2.0, interpolate=True)
    ptr_polygon_ordinate = visual.Line(
        win=win, name='ptr_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    ptr_polygon_petey_midpoint = visual.Line(
        win=win, name='ptr_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor='red', fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-4.0, interpolate=True)
    ptr_polygon_left_bound = visual.Line(
        win=win, name='ptr_polygon_left_bound',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor='white', fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-5.0, interpolate=True)
    ptr_polygon_right_bound = visual.Line(
        win=win, name='ptr_polygon_right_bound',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor='white', fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-6.0, interpolate=True)
    ptr_image_pigeon = visual.ImageStim(
        win=win,
        name='ptr_image_pigeon', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    ptr_image_right_seeds = visual.ImageStim(
        win=win,
        name='ptr_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-8.0)
    ptr_image_left_seeds = visual.ImageStim(
        win=win,
        name='ptr_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-9.0)
    ptr_polygon_coins_bar = visual.Rect(
        win=win, name='ptr_polygon_coins_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor='white', fillColor='white',
        opacity=None, depth=-10.0, interpolate=True)
    ptr_polygon_steps_bar = visual.Rect(
        win=win, name='ptr_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor='white', fillColor='white',
        opacity=None, depth=-11.0, interpolate=True)
    ptr_polygon_status_bar_apex = visual.Line(
        win=win, name='ptr_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor='black', fillColor='black',
        opacity=None, depth=-12.0, interpolate=True)
    ptr_polygon_status_bar_base = visual.Line(
        win=win, name='ptr_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-13.0, interpolate=True)
    ptr_text_coins = visual.TextStim(win=win, name='ptr_text_coins',
        text='',
        font='Open Sans',
        pos=(COINS_TEXT_X,COINS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-14.0);
    ptr_text_coins_label = visual.TextStim(win=win, name='ptr_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-15.0);
    ptr_text_coins_properties = visual.TextStim(win=win, name='ptr_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-16.0);
    ptr_text_steps = visual.TextStim(win=win, name='ptr_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-17.0);
    ptr_text_steps_label = visual.TextStim(win=win, name='ptr_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X,STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color='black', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-18.0);
    ptr_text_steps_properties = visual.TextStim(win=win, name='ptr_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-19.0);
    ptr_text_instruction = visual.TextStim(win=win, name='ptr_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-20.0);
    
    # --- Initialize components for Routine "show_message" ---
    sm_kbd_continue = keyboard.Keyboard(deviceName='sm_kbd_continue')
    # Run 'Begin Experiment' code from sm_code_update
    box_width = 0
    box_height = 0
    box_x = 0
    box_y = 0
    sm_polygon_box = visual.Rect(
        win=win, name='sm_polygon_box',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-2.0, interpolate=True)
    sm_polygon_abscissa = visual.Line(
        win=win, name='sm_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-3.0, interpolate=True)
    sm_polygon_ordinate = visual.Line(
        win=win, name='sm_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-4.0, interpolate=True)
    sm_polygon_petey_midpoint = visual.Line(
        win=win, name='sm_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-5.0, interpolate=True)
    sm_image_right_seeds = visual.ImageStim(
        win=win,
        name='sm_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-6.0)
    sm_image_left_seeds = visual.ImageStim(
        win=win,
        name='sm_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-7.0)
    sm_image_petey = visual.ImageStim(
        win=win,
        name='sm_image_petey', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-8.0)
    sm_polygon_coin_bar = visual.Rect(
        win=win, name='sm_polygon_coin_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-9.0, interpolate=True)
    sm_polygon_steps_bar = visual.Rect(
        win=win, name='sm_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-10.0, interpolate=True)
    sm_polygon_status_bar_apex = visual.Line(
        win=win, name='sm_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-11.0, interpolate=True)
    sm_polygon_status_bar_base = visual.Line(
        win=win, name='sm_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-12.0, interpolate=True)
    sm_text_coins = visual.TextStim(win=win, name='sm_text_coins',
        text='',
        font='Open Sans',
        pos=(COINS_TEXT_X,COINS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-13.0);
    sm_text_coins_label = visual.TextStim(win=win, name='sm_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-14.0);
    sm_text_coins_properties = visual.TextStim(win=win, name='sm_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-15.0);
    sm_text_steps = visual.TextStim(win=win, name='sm_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-16.0);
    sm_text_steps_label = visual.TextStim(win=win, name='sm_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-17.0);
    sm_text_steps_properties = visual.TextStim(win=win, name='sm_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-18.0);
    sm_text_message = visual.TextStim(win=win, name='sm_text_message',
        text='',
        font='Open Sans',
        pos=(MESSAGE_X,MESSAGE_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-19.0);
    sm_text_instruction = visual.TextStim(win=win, name='sm_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-20.0);
    
    # --- Initialize components for Routine "final_message" ---
    
    # --- Initialize components for Routine "show_final_message" ---
    sfm_kbd_continue = keyboard.Keyboard(deviceName='sfm_kbd_continue')
    # Run 'Begin Experiment' code from sfm_code_update
    box_width = 0
    box_height = 0
    box_x = 0
    box_y = 0
    sfm_text_message = visual.TextStim(win=win, name='sfm_text_message',
        text='',
        font='Open Sans',
        pos=(MESSAGE_X,MESSAGE_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-2.0);
    sfm_text_steps_properties = visual.TextStim(win=win, name='sfm_text_steps_properties',
        text='',
        font='Open Sans',
        pos=(STEPS_PROPERTIES_TEXT_X,STEPS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-3.0);
    sfm_text_steps_label = visual.TextStim(win=win, name='sfm_text_steps_label',
        text='Steps remaining',
        font='Open Sans',
        pos=(STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-4.0);
    sfm_text_steps = visual.TextStim(win=win, name='sfm_text_steps',
        text='',
        font='Open Sans',
        pos=(STEPS_TEXT_X,STEPS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-5.0);
    sfm_text_coins_properties = visual.TextStim(win=win, name='sfm_text_coins_properties',
        text='',
        font='Open Sans',
        pos=(COINS_PROPERTIES_TEXT_X,COINS_PROPERTIES_TEXT_Y), draggable=False, height=0.025, wrapWidth=None, ori=0.0, 
        color=[1.0000, 1.0000, 1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-6.0);
    sfm_text_coins_label = visual.TextStim(win=win, name='sfm_text_coins_label',
        text='Coins earned',
        font='Open Sans',
        pos=(COINS_LABEL_TEXT_X,COINS_LABEL_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-7.0);
    sfm_text_coins = visual.TextStim(win=win, name='sfm_text_coins',
        text='',
        font='Open Sans',
        pos=(COINS_TEXT_X,COINS_TEXT_Y), draggable=False, height=0.03, wrapWidth=None, ori=0.0, 
        color=[-1.0000, -1.0000, -1.0000], colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-8.0);
    sfm_polygon_status_bar_base = visual.Line(
        win=win, name='sfm_polygon_status_bar_base',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-9.0, interpolate=True)
    sfm_polygon_status_bar_apex = visual.Line(
        win=win, name='sfm_polygon_status_bar_apex',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-10.0, interpolate=True)
    sfm_polygon_steps_bar = visual.Rect(
        win=win, name='sfm_polygon_steps_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-11.0, interpolate=True)
    sfm_polygon_coin_bar = visual.Rect(
        win=win, name='sfm_polygon_coin_bar',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=1.0,
        colorSpace='rgb', lineColor=[1.0000, 1.0000, 1.0000], fillColor=[1.0000, 1.0000, 1.0000],
        opacity=None, depth=-12.0, interpolate=True)
    sfm_image_petey = visual.ImageStim(
        win=win,
        name='sfm_image_petey', 
        image='Resources/pigeon.png', mask=None, anchor='center',
        ori=0.0, pos=[0,0], draggable=False, size=1.0,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-13.0)
    sfm_image_left_seeds = visual.ImageStim(
        win=win,
        name='sfm_image_left_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(-EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-14.0)
    sfm_image_right_seeds = visual.ImageStim(
        win=win,
        name='sfm_image_right_seeds', 
        image='Resources/seeds.png', mask=None, anchor='center',
        ori=0.0, pos=(EDGE_DISTANCE,ORIGIN_Y), draggable=False, size=(SEEDS_WIDTH,SEEDS_HEIGHT),
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-15.0)
    sfm_polygon_petey_midpoint = visual.Line(
        win=win, name='sfm_polygon_petey_midpoint',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-16.0, interpolate=True)
    sfm_polygon_ordinate = visual.Line(
        win=win, name='sfm_polygon_ordinate',
        size=[1.0, 1.0],
        ori=90.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-17.0, interpolate=True)
    sfm_polygon_abscissa = visual.Line(
        win=win, name='sfm_polygon_abscissa',
        size=[1.0, 1.0],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor=[-1.0000, -1.0000, -1.0000], fillColor=[-1.0000, -1.0000, -1.0000],
        opacity=None, depth=-18.0, interpolate=True)
    sfm_polygon_box = visual.Rect(
        win=win, name='sfm_polygon_box',
        width=[1.0, 1.0][0], height=[1.0, 1.0][1],
        ori=0.0, pos=[0,0], draggable=False, anchor='center',
        lineWidth=3.0,
        colorSpace='rgb', lineColor=[1.0000, -1.0000, -1.0000], fillColor=[0.0000, 0.0000, 0.0000],
        opacity=None, depth=-19.0, interpolate=True)
    sfm_text_instruction = visual.TextStim(win=win, name='sfm_text_instruction',
        text='',
        font='Open Sans',
        pos=(INSTRUCTION_X,INSTRUCTION_Y), draggable=False, height=0.02, wrapWidth=None, ori=0.0, 
        color=INSTRUCTION_COLOR, colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-20.0);
    
    # create some handy timers
    
    # global clock to track the time since experiment started
    if globalClock is None:
        # create a clock if not given one
        globalClock = core.Clock()
    if isinstance(globalClock, str):
        # if given a string, make a clock accoridng to it
        if globalClock == 'float':
            # get timestamps as a simple value
            globalClock = core.Clock(format='float')
        elif globalClock == 'iso':
            # get timestamps in ISO format
            globalClock = core.Clock(format='%Y-%m-%d_%H:%M:%S.%f%z')
        else:
            # get timestamps in a custom format
            globalClock = core.Clock(format=globalClock)
    if ioServer is not None:
        ioServer.syncClock(globalClock)
    logging.setDefaultClock(globalClock)
    # routine timer to track time remaining of each (possibly non-slip) routine
    routineTimer = core.Clock()
    win.flip()  # flip window to reset last flip timer
    # store the exact time the global clock started
    expInfo['expStart'] = data.getDateStr(
        format='%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6
    )
    
    # --- Prepare to start Routine "consent" ---
    # create an object to store info about Routine consent
    consent = data.Routine(
        name='consent',
        components=[consent_txt, consent_box, consent_mouse, consent_acknowledge],
    )
    consent.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    consent_txt.reset()
    # setup some python lists for storing info about the consent_mouse
    consent_mouse.clicked_name = []
    gotValidClick = False  # until a click is received
    consent_acknowledge.reset()
    # store start times for consent
    consent.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    consent.tStart = globalClock.getTime(format='float')
    consent.status = STARTED
    thisExp.addData('consent.started', consent.tStart)
    consent.maxDuration = None
    # keep track of which components have finished
    consentComponents = consent.components
    for thisComponent in consent.components:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "consent" ---
    consent.forceEnded = routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *consent_txt* updates
        
        # if consent_txt is starting this frame...
        if consent_txt.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            consent_txt.frameNStart = frameN  # exact frame index
            consent_txt.tStart = t  # local t and not account for scr refresh
            consent_txt.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(consent_txt, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'consent_txt.started')
            # update status
            consent_txt.status = STARTED
            consent_txt.setAutoDraw(True)
        
        # if consent_txt is active this frame...
        if consent_txt.status == STARTED:
            # update params
            pass
        
        # *consent_box* updates
        
        # if consent_box is starting this frame...
        if consent_box.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            consent_box.frameNStart = frameN  # exact frame index
            consent_box.tStart = t  # local t and not account for scr refresh
            consent_box.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(consent_box, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'consent_box.started')
            # update status
            consent_box.status = STARTED
            consent_box.setAutoDraw(True)
        
        # if consent_box is active this frame...
        if consent_box.status == STARTED:
            # update params
            pass
        # *consent_mouse* updates
        
        # if consent_mouse is starting this frame...
        if consent_mouse.status == NOT_STARTED and t >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            consent_mouse.frameNStart = frameN  # exact frame index
            consent_mouse.tStart = t  # local t and not account for scr refresh
            consent_mouse.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(consent_mouse, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.addData('consent_mouse.started', t)
            # update status
            consent_mouse.status = STARTED
            consent_mouse.mouseClock.reset()
            prevButtonState = consent_mouse.getPressed()  # if button is down already this ISN'T a new click
        if consent_mouse.status == STARTED:  # only update if started and not finished!
            buttons = consent_mouse.getPressed()
            if buttons != prevButtonState:  # button state changed?
                prevButtonState = buttons
                if sum(buttons) > 0:  # state changed to a new click
                    # check if the mouse was inside our 'clickable' objects
                    gotValidClick = False
                    clickableList = environmenttools.getFromNames(consent_box, namespace=locals())
                    for obj in clickableList:
                        # is this object clicked on?
                        if obj.contains(consent_mouse):
                            gotValidClick = True
                            consent_mouse.clicked_name.append(obj.name)
                    if not gotValidClick:
                        consent_mouse.clicked_name.append(None)
                    if gotValidClick:  
                        continueRoutine = False  # end routine on response
        
        # *consent_acknowledge* updates
        
        # if consent_acknowledge is starting this frame...
        if consent_acknowledge.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            consent_acknowledge.frameNStart = frameN  # exact frame index
            consent_acknowledge.tStart = t  # local t and not account for scr refresh
            consent_acknowledge.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(consent_acknowledge, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'consent_acknowledge.started')
            # update status
            consent_acknowledge.status = STARTED
            consent_acknowledge.setAutoDraw(True)
        
        # if consent_acknowledge is active this frame...
        if consent_acknowledge.status == STARTED:
            # update params
            pass
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, win=win)
            return
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
            )
            # skip the frame we paused on
            continue
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            consent.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in consent.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "consent" ---
    for thisComponent in consent.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for consent
    consent.tStop = globalClock.getTime(format='float')
    consent.tStopRefresh = tThisFlipGlobal
    thisExp.addData('consent.stopped', consent.tStop)
    # store data for thisExp (ExperimentHandler)
    x, y = consent_mouse.getPos()
    buttons = consent_mouse.getPressed()
    if sum(buttons):
        # check if the mouse was inside our 'clickable' objects
        gotValidClick = False
        clickableList = environmenttools.getFromNames(consent_box, namespace=locals())
        for obj in clickableList:
            # is this object clicked on?
            if obj.contains(consent_mouse):
                gotValidClick = True
                consent_mouse.clicked_name.append(obj.name)
        if not gotValidClick:
            consent_mouse.clicked_name.append(None)
    thisExp.addData('consent_mouse.x', x)
    thisExp.addData('consent_mouse.y', y)
    thisExp.addData('consent_mouse.leftButton', buttons[0])
    thisExp.addData('consent_mouse.midButton', buttons[1])
    thisExp.addData('consent_mouse.rightButton', buttons[2])
    if len(consent_mouse.clicked_name):
        thisExp.addData('consent_mouse.clicked_name', consent_mouse.clicked_name[0])
    thisExp.nextEntry()
    # the Routine "consent" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "experiment_setup" ---
    # create an object to store info about Routine experiment_setup
    experiment_setup = data.Routine(
        name='experiment_setup',
        components=[],
    )
    experiment_setup.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # store start times for experiment_setup
    experiment_setup.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    experiment_setup.tStart = globalClock.getTime(format='float')
    experiment_setup.status = STARTED
    thisExp.addData('experiment_setup.started', experiment_setup.tStart)
    experiment_setup.maxDuration = None
    # keep track of which components have finished
    experiment_setupComponents = experiment_setup.components
    for thisComponent in experiment_setup.components:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "experiment_setup" ---
    experiment_setup.forceEnded = routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, win=win)
            return
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
            )
            # skip the frame we paused on
            continue
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            experiment_setup.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in experiment_setup.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "experiment_setup" ---
    for thisComponent in experiment_setup.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for experiment_setup
    experiment_setup.tStop = globalClock.getTime(format='float')
    experiment_setup.tStopRefresh = tThisFlipGlobal
    thisExp.addData('experiment_setup.stopped', experiment_setup.tStop)
    thisExp.nextEntry()
    # the Routine "experiment_setup" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    general_instructions_loop = data.TrialHandler2(
        name='general_instructions_loop',
        nReps=1.0, 
        method='sequential', 
        extraInfo=expInfo, 
        originPath=-1, 
        trialList=data.importConditions('Resources/PigeonGeneralInstructionsConditions.xlsx'), 
        seed=None, 
    )
    thisExp.addLoop(general_instructions_loop)  # add the loop to the experiment
    thisGeneral_instructions_loop = general_instructions_loop.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisGeneral_instructions_loop.rgb)
    if thisGeneral_instructions_loop != None:
        for paramName in thisGeneral_instructions_loop:
            globals()[paramName] = thisGeneral_instructions_loop[paramName]
    
    for thisGeneral_instructions_loop in general_instructions_loop:
        currentLoop = general_instructions_loop
        thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
        # abbreviate parameter names if possible (e.g. rgb = thisGeneral_instructions_loop.rgb)
        if thisGeneral_instructions_loop != None:
            for paramName in thisGeneral_instructions_loop:
                globals()[paramName] = thisGeneral_instructions_loop[paramName]
        
        # --- Prepare to start Routine "general_instructions_setup" ---
        # create an object to store info about Routine general_instructions_setup
        general_instructions_setup = data.Routine(
            name='general_instructions_setup',
            components=[],
        )
        general_instructions_setup.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # Run 'Begin Routine' code from gis_code_setup
        status_show = False
        box_show = False
        pigeon_show = gic_pigeon_show
        seeds_show = gic_seeds_show
        message_string = gic_message
        # store start times for general_instructions_setup
        general_instructions_setup.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        general_instructions_setup.tStart = globalClock.getTime(format='float')
        general_instructions_setup.status = STARTED
        thisExp.addData('general_instructions_setup.started', general_instructions_setup.tStart)
        general_instructions_setup.maxDuration = None
        # keep track of which components have finished
        general_instructions_setupComponents = general_instructions_setup.components
        for thisComponent in general_instructions_setup.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "general_instructions_setup" ---
        # if trial has changed, end Routine now
        if isinstance(general_instructions_loop, data.TrialHandler2) and thisGeneral_instructions_loop.thisN != general_instructions_loop.thisTrial.thisN:
            continueRoutine = False
        general_instructions_setup.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                general_instructions_setup.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in general_instructions_setup.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "general_instructions_setup" ---
        for thisComponent in general_instructions_setup.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for general_instructions_setup
        general_instructions_setup.tStop = globalClock.getTime(format='float')
        general_instructions_setup.tStopRefresh = tThisFlipGlobal
        thisExp.addData('general_instructions_setup.stopped', general_instructions_setup.tStop)
        # the Routine "general_instructions_setup" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "show_message" ---
        # create an object to store info about Routine show_message
        show_message = data.Routine(
            name='show_message',
            components=[sm_kbd_continue, sm_polygon_box, sm_polygon_abscissa, sm_polygon_ordinate, sm_polygon_petey_midpoint, sm_image_right_seeds, sm_image_left_seeds, sm_image_petey, sm_polygon_coin_bar, sm_polygon_steps_bar, sm_polygon_status_bar_apex, sm_polygon_status_bar_base, sm_text_coins, sm_text_coins_label, sm_text_coins_properties, sm_text_steps, sm_text_steps_label, sm_text_steps_properties, sm_text_message, sm_text_instruction],
        )
        show_message.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # create starting attributes for sm_kbd_continue
        sm_kbd_continue.keys = []
        sm_kbd_continue.rt = []
        _sm_kbd_continue_allKeys = []
        # Run 'Begin Routine' code from sm_code_update
        # Clear all keyboard events
        sm_kbd_continue.clearEvents()
        sm_polygon_box.setPos((box_x,box_y))
        sm_polygon_box.setSize((box_width,box_height))
        sm_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
        sm_polygon_abscissa.setSize((ABSCISSA_WIDTH,LINE_HEIGHT))
        sm_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
        sm_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
        sm_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
        sm_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
        sm_image_petey.setPos((pigeon_shown_x,ORIGIN_Y))
        sm_image_petey.setSize((pigeon_flip, 0.1))
        sm_polygon_coin_bar.setPos((COINS_BAR_X,coins_bar_center))
        sm_polygon_coin_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
        sm_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
        sm_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
        sm_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
        sm_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
        sm_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
        sm_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
        sm_text_coins.setText(coins_count)
        sm_text_coins_properties.setText(coins_string)
        sm_text_steps.setText(steps_max - steps_count)
        sm_text_steps_properties.setText(steps_string)
        sm_text_message.setText(message_string)
        sm_text_instruction.setText(CONTINUE_STRING)
        # store start times for show_message
        show_message.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        show_message.tStart = globalClock.getTime(format='float')
        show_message.status = STARTED
        thisExp.addData('show_message.started', show_message.tStart)
        show_message.maxDuration = None
        # keep track of which components have finished
        show_messageComponents = show_message.components
        for thisComponent in show_message.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "show_message" ---
        # if trial has changed, end Routine now
        if isinstance(general_instructions_loop, data.TrialHandler2) and thisGeneral_instructions_loop.thisN != general_instructions_loop.thisTrial.thisN:
            continueRoutine = False
        show_message.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *sm_kbd_continue* updates
            waitOnFlip = False
            
            # if sm_kbd_continue is starting this frame...
            if sm_kbd_continue.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sm_kbd_continue.frameNStart = frameN  # exact frame index
                sm_kbd_continue.tStart = t  # local t and not account for scr refresh
                sm_kbd_continue.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_kbd_continue, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_kbd_continue.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(sm_kbd_continue.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(sm_kbd_continue.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if sm_kbd_continue.status == STARTED and not waitOnFlip:
                theseKeys = sm_kbd_continue.getKeys(keyList=['space', 'up', 'down', 'left', 'right'], ignoreKeys=["escape"], waitRelease=False)
                _sm_kbd_continue_allKeys.extend(theseKeys)
                if len(_sm_kbd_continue_allKeys):
                    sm_kbd_continue.keys = _sm_kbd_continue_allKeys[-1].name  # just the last key pressed
                    sm_kbd_continue.rt = _sm_kbd_continue_allKeys[-1].rt
                    sm_kbd_continue.duration = _sm_kbd_continue_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            
            # *sm_polygon_box* updates
            
            # if sm_polygon_box is starting this frame...
            if sm_polygon_box.status == NOT_STARTED and box_show:
                # keep track of start time/frame for later
                sm_polygon_box.frameNStart = frameN  # exact frame index
                sm_polygon_box.tStart = t  # local t and not account for scr refresh
                sm_polygon_box.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_box, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_box.status = STARTED
                sm_polygon_box.setAutoDraw(True)
            
            # if sm_polygon_box is active this frame...
            if sm_polygon_box.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_abscissa* updates
            
            # if sm_polygon_abscissa is starting this frame...
            if sm_polygon_abscissa.status == NOT_STARTED and pigeon_show:
                # keep track of start time/frame for later
                sm_polygon_abscissa.frameNStart = frameN  # exact frame index
                sm_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                sm_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_abscissa.status = STARTED
                sm_polygon_abscissa.setAutoDraw(True)
            
            # if sm_polygon_abscissa is active this frame...
            if sm_polygon_abscissa.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_ordinate* updates
            
            # if sm_polygon_ordinate is starting this frame...
            if sm_polygon_ordinate.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_ordinate.frameNStart = frameN  # exact frame index
                sm_polygon_ordinate.tStart = t  # local t and not account for scr refresh
                sm_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_ordinate.status = STARTED
                sm_polygon_ordinate.setAutoDraw(True)
            
            # if sm_polygon_ordinate is active this frame...
            if sm_polygon_ordinate.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_petey_midpoint* updates
            
            # if sm_polygon_petey_midpoint is starting this frame...
            if sm_polygon_petey_midpoint.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
                sm_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
                sm_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_petey_midpoint.status = STARTED
                sm_polygon_petey_midpoint.setAutoDraw(True)
            
            # if sm_polygon_petey_midpoint is active this frame...
            if sm_polygon_petey_midpoint.status == STARTED:
                # update params
                pass
            
            # *sm_image_right_seeds* updates
            
            # if sm_image_right_seeds is starting this frame...
            if sm_image_right_seeds.status == NOT_STARTED and seeds_show:
                # keep track of start time/frame for later
                sm_image_right_seeds.frameNStart = frameN  # exact frame index
                sm_image_right_seeds.tStart = t  # local t and not account for scr refresh
                sm_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_right_seeds.status = STARTED
                sm_image_right_seeds.setAutoDraw(True)
            
            # if sm_image_right_seeds is active this frame...
            if sm_image_right_seeds.status == STARTED:
                # update params
                pass
            
            # *sm_image_left_seeds* updates
            
            # if sm_image_left_seeds is starting this frame...
            if sm_image_left_seeds.status == NOT_STARTED and seeds_show:
                # keep track of start time/frame for later
                sm_image_left_seeds.frameNStart = frameN  # exact frame index
                sm_image_left_seeds.tStart = t  # local t and not account for scr refresh
                sm_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_left_seeds.status = STARTED
                sm_image_left_seeds.setAutoDraw(True)
            
            # if sm_image_left_seeds is active this frame...
            if sm_image_left_seeds.status == STARTED:
                # update params
                pass
            
            # *sm_image_petey* updates
            
            # if sm_image_petey is starting this frame...
            if sm_image_petey.status == NOT_STARTED and pigeon_show:
                # keep track of start time/frame for later
                sm_image_petey.frameNStart = frameN  # exact frame index
                sm_image_petey.tStart = t  # local t and not account for scr refresh
                sm_image_petey.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_petey, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_petey.status = STARTED
                sm_image_petey.setAutoDraw(True)
            
            # if sm_image_petey is active this frame...
            if sm_image_petey.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_coin_bar* updates
            
            # if sm_polygon_coin_bar is starting this frame...
            if sm_polygon_coin_bar.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_coin_bar.frameNStart = frameN  # exact frame index
                sm_polygon_coin_bar.tStart = t  # local t and not account for scr refresh
                sm_polygon_coin_bar.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_coin_bar, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_coin_bar.status = STARTED
                sm_polygon_coin_bar.setAutoDraw(True)
            
            # if sm_polygon_coin_bar is active this frame...
            if sm_polygon_coin_bar.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_steps_bar* updates
            
            # if sm_polygon_steps_bar is starting this frame...
            if sm_polygon_steps_bar.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_steps_bar.frameNStart = frameN  # exact frame index
                sm_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
                sm_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_steps_bar.status = STARTED
                sm_polygon_steps_bar.setAutoDraw(True)
            
            # if sm_polygon_steps_bar is active this frame...
            if sm_polygon_steps_bar.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_status_bar_apex* updates
            
            # if sm_polygon_status_bar_apex is starting this frame...
            if sm_polygon_status_bar_apex.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
                sm_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
                sm_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_status_bar_apex.status = STARTED
                sm_polygon_status_bar_apex.setAutoDraw(True)
            
            # if sm_polygon_status_bar_apex is active this frame...
            if sm_polygon_status_bar_apex.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_status_bar_base* updates
            
            # if sm_polygon_status_bar_base is starting this frame...
            if sm_polygon_status_bar_base.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_status_bar_base.frameNStart = frameN  # exact frame index
                sm_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
                sm_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_status_bar_base.status = STARTED
                sm_polygon_status_bar_base.setAutoDraw(True)
            
            # if sm_polygon_status_bar_base is active this frame...
            if sm_polygon_status_bar_base.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins* updates
            
            # if sm_text_coins is starting this frame...
            if sm_text_coins.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins.frameNStart = frameN  # exact frame index
                sm_text_coins.tStart = t  # local t and not account for scr refresh
                sm_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins.status = STARTED
                sm_text_coins.setAutoDraw(True)
            
            # if sm_text_coins is active this frame...
            if sm_text_coins.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins_label* updates
            
            # if sm_text_coins_label is starting this frame...
            if sm_text_coins_label.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins_label.frameNStart = frameN  # exact frame index
                sm_text_coins_label.tStart = t  # local t and not account for scr refresh
                sm_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins_label, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins_label.status = STARTED
                sm_text_coins_label.setAutoDraw(True)
            
            # if sm_text_coins_label is active this frame...
            if sm_text_coins_label.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins_properties* updates
            
            # if sm_text_coins_properties is starting this frame...
            if sm_text_coins_properties.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins_properties.frameNStart = frameN  # exact frame index
                sm_text_coins_properties.tStart = t  # local t and not account for scr refresh
                sm_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins_properties.status = STARTED
                sm_text_coins_properties.setAutoDraw(True)
            
            # if sm_text_coins_properties is active this frame...
            if sm_text_coins_properties.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps* updates
            
            # if sm_text_steps is starting this frame...
            if sm_text_steps.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps.frameNStart = frameN  # exact frame index
                sm_text_steps.tStart = t  # local t and not account for scr refresh
                sm_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps.status = STARTED
                sm_text_steps.setAutoDraw(True)
            
            # if sm_text_steps is active this frame...
            if sm_text_steps.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps_label* updates
            
            # if sm_text_steps_label is starting this frame...
            if sm_text_steps_label.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps_label.frameNStart = frameN  # exact frame index
                sm_text_steps_label.tStart = t  # local t and not account for scr refresh
                sm_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps_label, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps_label.status = STARTED
                sm_text_steps_label.setAutoDraw(True)
            
            # if sm_text_steps_label is active this frame...
            if sm_text_steps_label.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps_properties* updates
            
            # if sm_text_steps_properties is starting this frame...
            if sm_text_steps_properties.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps_properties.frameNStart = frameN  # exact frame index
                sm_text_steps_properties.tStart = t  # local t and not account for scr refresh
                sm_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps_properties.status = STARTED
                sm_text_steps_properties.setAutoDraw(True)
            
            # if sm_text_steps_properties is active this frame...
            if sm_text_steps_properties.status == STARTED:
                # update params
                pass
            
            # *sm_text_message* updates
            
            # if sm_text_message is starting this frame...
            if sm_text_message.status == NOT_STARTED and message_show:
                # keep track of start time/frame for later
                sm_text_message.frameNStart = frameN  # exact frame index
                sm_text_message.tStart = t  # local t and not account for scr refresh
                sm_text_message.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_message, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_message.status = STARTED
                sm_text_message.setAutoDraw(True)
            
            # if sm_text_message is active this frame...
            if sm_text_message.status == STARTED:
                # update params
                pass
            
            # *sm_text_instruction* updates
            
            # if sm_text_instruction is starting this frame...
            if sm_text_instruction.status == NOT_STARTED and instruction_show:
                # keep track of start time/frame for later
                sm_text_instruction.frameNStart = frameN  # exact frame index
                sm_text_instruction.tStart = t  # local t and not account for scr refresh
                sm_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_instruction, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_instruction.status = STARTED
                sm_text_instruction.setAutoDraw(True)
            
            # if sm_text_instruction is active this frame...
            if sm_text_instruction.status == STARTED:
                # update params
                pass
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                show_message.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in show_message.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "show_message" ---
        for thisComponent in show_message.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for show_message
        show_message.tStop = globalClock.getTime(format='float')
        show_message.tStopRefresh = tThisFlipGlobal
        thisExp.addData('show_message.stopped', show_message.tStop)
        # Run 'End Routine' code from sm_code_update
        # Ick, but whatever. Makes sure that we don't ever
        # drop unexpectedly out of a move loop without 
        # resetting
        trial_start = True
        # the Routine "show_message" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
    # completed 1.0 repeats of 'general_instructions_loop'
    
    
    # set up handler to look after randomisation of conditions etc
    simple_demo_block_loop = data.TrialHandler2(
        name='simple_demo_block_loop',
        nReps=1.0, 
        method='sequential', 
        extraInfo=expInfo, 
        originPath=-1, 
        trialList=data.importConditions('Resources/PigeonSimpleDemoConditions.xlsx'), 
        seed=None, 
    )
    thisExp.addLoop(simple_demo_block_loop)  # add the loop to the experiment
    thisSimple_demo_block_loop = simple_demo_block_loop.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisSimple_demo_block_loop.rgb)
    if thisSimple_demo_block_loop != None:
        for paramName in thisSimple_demo_block_loop:
            globals()[paramName] = thisSimple_demo_block_loop[paramName]
    
    for thisSimple_demo_block_loop in simple_demo_block_loop:
        currentLoop = simple_demo_block_loop
        thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
        # abbreviate parameter names if possible (e.g. rgb = thisSimple_demo_block_loop.rgb)
        if thisSimple_demo_block_loop != None:
            for paramName in thisSimple_demo_block_loop:
                globals()[paramName] = thisSimple_demo_block_loop[paramName]
        
        # --- Prepare to start Routine "simple_demo_setup" ---
        # create an object to store info about Routine simple_demo_setup
        simple_demo_setup = data.Routine(
            name='simple_demo_setup',
            components=[],
        )
        simple_demo_setup.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # Run 'Begin Routine' code from sds_code_setup
        # Set block conditions
        seeds_show = True
        pigeon_show = True
        status_show = False
        box_show = False
        demo_start = True
        
        # Get variables from conditions file
        message_string = sdc_message
        step_mean = sdc_step_mean
        step_std = sdc_step_std
        # store start times for simple_demo_setup
        simple_demo_setup.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        simple_demo_setup.tStart = globalClock.getTime(format='float')
        simple_demo_setup.status = STARTED
        thisExp.addData('simple_demo_setup.started', simple_demo_setup.tStart)
        simple_demo_setup.maxDuration = None
        # keep track of which components have finished
        simple_demo_setupComponents = simple_demo_setup.components
        for thisComponent in simple_demo_setup.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "simple_demo_setup" ---
        # if trial has changed, end Routine now
        if isinstance(simple_demo_block_loop, data.TrialHandler2) and thisSimple_demo_block_loop.thisN != simple_demo_block_loop.thisTrial.thisN:
            continueRoutine = False
        simple_demo_setup.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                simple_demo_setup.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in simple_demo_setup.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "simple_demo_setup" ---
        for thisComponent in simple_demo_setup.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for simple_demo_setup
        simple_demo_setup.tStop = globalClock.getTime(format='float')
        simple_demo_setup.tStopRefresh = tThisFlipGlobal
        thisExp.addData('simple_demo_setup.stopped', simple_demo_setup.tStop)
        # the Routine "simple_demo_setup" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "show_message" ---
        # create an object to store info about Routine show_message
        show_message = data.Routine(
            name='show_message',
            components=[sm_kbd_continue, sm_polygon_box, sm_polygon_abscissa, sm_polygon_ordinate, sm_polygon_petey_midpoint, sm_image_right_seeds, sm_image_left_seeds, sm_image_petey, sm_polygon_coin_bar, sm_polygon_steps_bar, sm_polygon_status_bar_apex, sm_polygon_status_bar_base, sm_text_coins, sm_text_coins_label, sm_text_coins_properties, sm_text_steps, sm_text_steps_label, sm_text_steps_properties, sm_text_message, sm_text_instruction],
        )
        show_message.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # create starting attributes for sm_kbd_continue
        sm_kbd_continue.keys = []
        sm_kbd_continue.rt = []
        _sm_kbd_continue_allKeys = []
        # Run 'Begin Routine' code from sm_code_update
        # Clear all keyboard events
        sm_kbd_continue.clearEvents()
        sm_polygon_box.setPos((box_x,box_y))
        sm_polygon_box.setSize((box_width,box_height))
        sm_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
        sm_polygon_abscissa.setSize((ABSCISSA_WIDTH,LINE_HEIGHT))
        sm_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
        sm_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
        sm_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
        sm_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
        sm_image_petey.setPos((pigeon_shown_x,ORIGIN_Y))
        sm_image_petey.setSize((pigeon_flip, 0.1))
        sm_polygon_coin_bar.setPos((COINS_BAR_X,coins_bar_center))
        sm_polygon_coin_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
        sm_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
        sm_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
        sm_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
        sm_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
        sm_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
        sm_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
        sm_text_coins.setText(coins_count)
        sm_text_coins_properties.setText(coins_string)
        sm_text_steps.setText(steps_max - steps_count)
        sm_text_steps_properties.setText(steps_string)
        sm_text_message.setText(message_string)
        sm_text_instruction.setText(CONTINUE_STRING)
        # store start times for show_message
        show_message.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        show_message.tStart = globalClock.getTime(format='float')
        show_message.status = STARTED
        thisExp.addData('show_message.started', show_message.tStart)
        show_message.maxDuration = None
        # keep track of which components have finished
        show_messageComponents = show_message.components
        for thisComponent in show_message.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "show_message" ---
        # if trial has changed, end Routine now
        if isinstance(simple_demo_block_loop, data.TrialHandler2) and thisSimple_demo_block_loop.thisN != simple_demo_block_loop.thisTrial.thisN:
            continueRoutine = False
        show_message.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *sm_kbd_continue* updates
            waitOnFlip = False
            
            # if sm_kbd_continue is starting this frame...
            if sm_kbd_continue.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sm_kbd_continue.frameNStart = frameN  # exact frame index
                sm_kbd_continue.tStart = t  # local t and not account for scr refresh
                sm_kbd_continue.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_kbd_continue, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_kbd_continue.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(sm_kbd_continue.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(sm_kbd_continue.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if sm_kbd_continue.status == STARTED and not waitOnFlip:
                theseKeys = sm_kbd_continue.getKeys(keyList=['space', 'up', 'down', 'left', 'right'], ignoreKeys=["escape"], waitRelease=False)
                _sm_kbd_continue_allKeys.extend(theseKeys)
                if len(_sm_kbd_continue_allKeys):
                    sm_kbd_continue.keys = _sm_kbd_continue_allKeys[-1].name  # just the last key pressed
                    sm_kbd_continue.rt = _sm_kbd_continue_allKeys[-1].rt
                    sm_kbd_continue.duration = _sm_kbd_continue_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            
            # *sm_polygon_box* updates
            
            # if sm_polygon_box is starting this frame...
            if sm_polygon_box.status == NOT_STARTED and box_show:
                # keep track of start time/frame for later
                sm_polygon_box.frameNStart = frameN  # exact frame index
                sm_polygon_box.tStart = t  # local t and not account for scr refresh
                sm_polygon_box.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_box, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_box.status = STARTED
                sm_polygon_box.setAutoDraw(True)
            
            # if sm_polygon_box is active this frame...
            if sm_polygon_box.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_abscissa* updates
            
            # if sm_polygon_abscissa is starting this frame...
            if sm_polygon_abscissa.status == NOT_STARTED and pigeon_show:
                # keep track of start time/frame for later
                sm_polygon_abscissa.frameNStart = frameN  # exact frame index
                sm_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                sm_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_abscissa.status = STARTED
                sm_polygon_abscissa.setAutoDraw(True)
            
            # if sm_polygon_abscissa is active this frame...
            if sm_polygon_abscissa.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_ordinate* updates
            
            # if sm_polygon_ordinate is starting this frame...
            if sm_polygon_ordinate.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_ordinate.frameNStart = frameN  # exact frame index
                sm_polygon_ordinate.tStart = t  # local t and not account for scr refresh
                sm_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_ordinate.status = STARTED
                sm_polygon_ordinate.setAutoDraw(True)
            
            # if sm_polygon_ordinate is active this frame...
            if sm_polygon_ordinate.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_petey_midpoint* updates
            
            # if sm_polygon_petey_midpoint is starting this frame...
            if sm_polygon_petey_midpoint.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
                sm_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
                sm_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_petey_midpoint.status = STARTED
                sm_polygon_petey_midpoint.setAutoDraw(True)
            
            # if sm_polygon_petey_midpoint is active this frame...
            if sm_polygon_petey_midpoint.status == STARTED:
                # update params
                pass
            
            # *sm_image_right_seeds* updates
            
            # if sm_image_right_seeds is starting this frame...
            if sm_image_right_seeds.status == NOT_STARTED and seeds_show:
                # keep track of start time/frame for later
                sm_image_right_seeds.frameNStart = frameN  # exact frame index
                sm_image_right_seeds.tStart = t  # local t and not account for scr refresh
                sm_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_right_seeds.status = STARTED
                sm_image_right_seeds.setAutoDraw(True)
            
            # if sm_image_right_seeds is active this frame...
            if sm_image_right_seeds.status == STARTED:
                # update params
                pass
            
            # *sm_image_left_seeds* updates
            
            # if sm_image_left_seeds is starting this frame...
            if sm_image_left_seeds.status == NOT_STARTED and seeds_show:
                # keep track of start time/frame for later
                sm_image_left_seeds.frameNStart = frameN  # exact frame index
                sm_image_left_seeds.tStart = t  # local t and not account for scr refresh
                sm_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_left_seeds.status = STARTED
                sm_image_left_seeds.setAutoDraw(True)
            
            # if sm_image_left_seeds is active this frame...
            if sm_image_left_seeds.status == STARTED:
                # update params
                pass
            
            # *sm_image_petey* updates
            
            # if sm_image_petey is starting this frame...
            if sm_image_petey.status == NOT_STARTED and pigeon_show:
                # keep track of start time/frame for later
                sm_image_petey.frameNStart = frameN  # exact frame index
                sm_image_petey.tStart = t  # local t and not account for scr refresh
                sm_image_petey.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_petey, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_petey.status = STARTED
                sm_image_petey.setAutoDraw(True)
            
            # if sm_image_petey is active this frame...
            if sm_image_petey.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_coin_bar* updates
            
            # if sm_polygon_coin_bar is starting this frame...
            if sm_polygon_coin_bar.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_coin_bar.frameNStart = frameN  # exact frame index
                sm_polygon_coin_bar.tStart = t  # local t and not account for scr refresh
                sm_polygon_coin_bar.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_coin_bar, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_coin_bar.status = STARTED
                sm_polygon_coin_bar.setAutoDraw(True)
            
            # if sm_polygon_coin_bar is active this frame...
            if sm_polygon_coin_bar.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_steps_bar* updates
            
            # if sm_polygon_steps_bar is starting this frame...
            if sm_polygon_steps_bar.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_steps_bar.frameNStart = frameN  # exact frame index
                sm_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
                sm_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_steps_bar.status = STARTED
                sm_polygon_steps_bar.setAutoDraw(True)
            
            # if sm_polygon_steps_bar is active this frame...
            if sm_polygon_steps_bar.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_status_bar_apex* updates
            
            # if sm_polygon_status_bar_apex is starting this frame...
            if sm_polygon_status_bar_apex.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
                sm_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
                sm_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_status_bar_apex.status = STARTED
                sm_polygon_status_bar_apex.setAutoDraw(True)
            
            # if sm_polygon_status_bar_apex is active this frame...
            if sm_polygon_status_bar_apex.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_status_bar_base* updates
            
            # if sm_polygon_status_bar_base is starting this frame...
            if sm_polygon_status_bar_base.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_status_bar_base.frameNStart = frameN  # exact frame index
                sm_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
                sm_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_status_bar_base.status = STARTED
                sm_polygon_status_bar_base.setAutoDraw(True)
            
            # if sm_polygon_status_bar_base is active this frame...
            if sm_polygon_status_bar_base.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins* updates
            
            # if sm_text_coins is starting this frame...
            if sm_text_coins.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins.frameNStart = frameN  # exact frame index
                sm_text_coins.tStart = t  # local t and not account for scr refresh
                sm_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins.status = STARTED
                sm_text_coins.setAutoDraw(True)
            
            # if sm_text_coins is active this frame...
            if sm_text_coins.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins_label* updates
            
            # if sm_text_coins_label is starting this frame...
            if sm_text_coins_label.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins_label.frameNStart = frameN  # exact frame index
                sm_text_coins_label.tStart = t  # local t and not account for scr refresh
                sm_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins_label, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins_label.status = STARTED
                sm_text_coins_label.setAutoDraw(True)
            
            # if sm_text_coins_label is active this frame...
            if sm_text_coins_label.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins_properties* updates
            
            # if sm_text_coins_properties is starting this frame...
            if sm_text_coins_properties.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins_properties.frameNStart = frameN  # exact frame index
                sm_text_coins_properties.tStart = t  # local t and not account for scr refresh
                sm_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins_properties.status = STARTED
                sm_text_coins_properties.setAutoDraw(True)
            
            # if sm_text_coins_properties is active this frame...
            if sm_text_coins_properties.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps* updates
            
            # if sm_text_steps is starting this frame...
            if sm_text_steps.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps.frameNStart = frameN  # exact frame index
                sm_text_steps.tStart = t  # local t and not account for scr refresh
                sm_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps.status = STARTED
                sm_text_steps.setAutoDraw(True)
            
            # if sm_text_steps is active this frame...
            if sm_text_steps.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps_label* updates
            
            # if sm_text_steps_label is starting this frame...
            if sm_text_steps_label.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps_label.frameNStart = frameN  # exact frame index
                sm_text_steps_label.tStart = t  # local t and not account for scr refresh
                sm_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps_label, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps_label.status = STARTED
                sm_text_steps_label.setAutoDraw(True)
            
            # if sm_text_steps_label is active this frame...
            if sm_text_steps_label.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps_properties* updates
            
            # if sm_text_steps_properties is starting this frame...
            if sm_text_steps_properties.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps_properties.frameNStart = frameN  # exact frame index
                sm_text_steps_properties.tStart = t  # local t and not account for scr refresh
                sm_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps_properties.status = STARTED
                sm_text_steps_properties.setAutoDraw(True)
            
            # if sm_text_steps_properties is active this frame...
            if sm_text_steps_properties.status == STARTED:
                # update params
                pass
            
            # *sm_text_message* updates
            
            # if sm_text_message is starting this frame...
            if sm_text_message.status == NOT_STARTED and message_show:
                # keep track of start time/frame for later
                sm_text_message.frameNStart = frameN  # exact frame index
                sm_text_message.tStart = t  # local t and not account for scr refresh
                sm_text_message.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_message, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_message.status = STARTED
                sm_text_message.setAutoDraw(True)
            
            # if sm_text_message is active this frame...
            if sm_text_message.status == STARTED:
                # update params
                pass
            
            # *sm_text_instruction* updates
            
            # if sm_text_instruction is starting this frame...
            if sm_text_instruction.status == NOT_STARTED and instruction_show:
                # keep track of start time/frame for later
                sm_text_instruction.frameNStart = frameN  # exact frame index
                sm_text_instruction.tStart = t  # local t and not account for scr refresh
                sm_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_instruction, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_instruction.status = STARTED
                sm_text_instruction.setAutoDraw(True)
            
            # if sm_text_instruction is active this frame...
            if sm_text_instruction.status == STARTED:
                # update params
                pass
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                show_message.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in show_message.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "show_message" ---
        for thisComponent in show_message.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for show_message
        show_message.tStop = globalClock.getTime(format='float')
        show_message.tStopRefresh = tThisFlipGlobal
        thisExp.addData('show_message.stopped', show_message.tStop)
        # Run 'End Routine' code from sm_code_update
        # Ick, but whatever. Makes sure that we don't ever
        # drop unexpectedly out of a move loop without 
        # resetting
        trial_start = True
        # the Routine "show_message" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # set up handler to look after randomisation of conditions etc
        simple_demo_trial_loop = data.TrialHandler2(
            name='simple_demo_trial_loop',
            nReps=1000.0, 
            method='sequential', 
            extraInfo=expInfo, 
            originPath=-1, 
            trialList=[None], 
            seed=None, 
        )
        thisExp.addLoop(simple_demo_trial_loop)  # add the loop to the experiment
        thisSimple_demo_trial_loop = simple_demo_trial_loop.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb = thisSimple_demo_trial_loop.rgb)
        if thisSimple_demo_trial_loop != None:
            for paramName in thisSimple_demo_trial_loop:
                globals()[paramName] = thisSimple_demo_trial_loop[paramName]
        
        for thisSimple_demo_trial_loop in simple_demo_trial_loop:
            currentLoop = simple_demo_trial_loop
            thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
            # abbreviate parameter names if possible (e.g. rgb = thisSimple_demo_trial_loop.rgb)
            if thisSimple_demo_trial_loop != None:
                for paramName in thisSimple_demo_trial_loop:
                    globals()[paramName] = thisSimple_demo_trial_loop[paramName]
            
            # --- Prepare to start Routine "simple_demo" ---
            # create an object to store info about Routine simple_demo
            simple_demo = data.Routine(
                name='simple_demo',
                components=[sd_polygon_abscissa, sd_image_left_seeds, sd_image_right_seeds, sd_image_pigeon],
            )
            simple_demo.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            # Run 'Begin Routine' code from sd_code_update
            if demo_start:
                # initialize pigeon variables
                pigeon_flip = 0.1
                pigeon_true_x = 0
                pigeon_shown_x = pigeon_true_x
                
                # Randomize direction, set std
                if randint(0,2)==0:
                    step_mean_val = -1*abs(step_mean)
                else:
                    step_mean_val = abs(step_mean)
                step_std_val = step_mean
                
                # unset flag
                demo_start = False
            else:
                # Check for completion
                if pigeon_true_x > -EDGE_DISTANCE and pigeon_true_x < EDGE_DISTANCE:
                    move_pigeon(step_mean_val, step_std_val)
                else:
                    simple_demo_trial_loop.finished = True
            sd_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
            sd_polygon_abscissa.setSize((EDGE_DISTANCE*2,LINE_HEIGHT))
            sd_image_pigeon.setPos((pigeon_shown_x,ORIGIN_Y))
            sd_image_pigeon.setSize((pigeon_flip, 0.1))
            # store start times for simple_demo
            simple_demo.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            simple_demo.tStart = globalClock.getTime(format='float')
            simple_demo.status = STARTED
            thisExp.addData('simple_demo.started', simple_demo.tStart)
            simple_demo.maxDuration = None
            # keep track of which components have finished
            simple_demoComponents = simple_demo.components
            for thisComponent in simple_demo.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "simple_demo" ---
            # if trial has changed, end Routine now
            if isinstance(simple_demo_trial_loop, data.TrialHandler2) and thisSimple_demo_trial_loop.thisN != simple_demo_trial_loop.thisTrial.thisN:
                continueRoutine = False
            simple_demo.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                # *sd_polygon_abscissa* updates
                
                # if sd_polygon_abscissa is starting this frame...
                if sd_polygon_abscissa.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    sd_polygon_abscissa.frameNStart = frameN  # exact frame index
                    sd_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                    sd_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sd_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sd_polygon_abscissa.status = STARTED
                    sd_polygon_abscissa.setAutoDraw(True)
                
                # if sd_polygon_abscissa is active this frame...
                if sd_polygon_abscissa.status == STARTED:
                    # update params
                    pass
                
                # if sd_polygon_abscissa is stopping this frame...
                if sd_polygon_abscissa.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > sd_polygon_abscissa.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                        # keep track of stop time/frame for later
                        sd_polygon_abscissa.tStop = t  # not accounting for scr refresh
                        sd_polygon_abscissa.tStopRefresh = tThisFlipGlobal  # on global time
                        sd_polygon_abscissa.frameNStop = frameN  # exact frame index
                        # update status
                        sd_polygon_abscissa.status = FINISHED
                        sd_polygon_abscissa.setAutoDraw(False)
                
                # *sd_image_left_seeds* updates
                
                # if sd_image_left_seeds is starting this frame...
                if sd_image_left_seeds.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    sd_image_left_seeds.frameNStart = frameN  # exact frame index
                    sd_image_left_seeds.tStart = t  # local t and not account for scr refresh
                    sd_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sd_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sd_image_left_seeds.status = STARTED
                    sd_image_left_seeds.setAutoDraw(True)
                
                # if sd_image_left_seeds is active this frame...
                if sd_image_left_seeds.status == STARTED:
                    # update params
                    pass
                
                # if sd_image_left_seeds is stopping this frame...
                if sd_image_left_seeds.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > sd_image_left_seeds.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                        # keep track of stop time/frame for later
                        sd_image_left_seeds.tStop = t  # not accounting for scr refresh
                        sd_image_left_seeds.tStopRefresh = tThisFlipGlobal  # on global time
                        sd_image_left_seeds.frameNStop = frameN  # exact frame index
                        # update status
                        sd_image_left_seeds.status = FINISHED
                        sd_image_left_seeds.setAutoDraw(False)
                
                # *sd_image_right_seeds* updates
                
                # if sd_image_right_seeds is starting this frame...
                if sd_image_right_seeds.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    sd_image_right_seeds.frameNStart = frameN  # exact frame index
                    sd_image_right_seeds.tStart = t  # local t and not account for scr refresh
                    sd_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sd_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sd_image_right_seeds.status = STARTED
                    sd_image_right_seeds.setAutoDraw(True)
                
                # if sd_image_right_seeds is active this frame...
                if sd_image_right_seeds.status == STARTED:
                    # update params
                    pass
                
                # if sd_image_right_seeds is stopping this frame...
                if sd_image_right_seeds.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > sd_image_right_seeds.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                        # keep track of stop time/frame for later
                        sd_image_right_seeds.tStop = t  # not accounting for scr refresh
                        sd_image_right_seeds.tStopRefresh = tThisFlipGlobal  # on global time
                        sd_image_right_seeds.frameNStop = frameN  # exact frame index
                        # update status
                        sd_image_right_seeds.status = FINISHED
                        sd_image_right_seeds.setAutoDraw(False)
                
                # *sd_image_pigeon* updates
                
                # if sd_image_pigeon is starting this frame...
                if sd_image_pigeon.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    sd_image_pigeon.frameNStart = frameN  # exact frame index
                    sd_image_pigeon.tStart = t  # local t and not account for scr refresh
                    sd_image_pigeon.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sd_image_pigeon, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sd_image_pigeon.status = STARTED
                    sd_image_pigeon.setAutoDraw(True)
                
                # if sd_image_pigeon is active this frame...
                if sd_image_pigeon.status == STARTED:
                    # update params
                    pass
                
                # if sd_image_pigeon is stopping this frame...
                if sd_image_pigeon.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > sd_image_pigeon.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                        # keep track of stop time/frame for later
                        sd_image_pigeon.tStop = t  # not accounting for scr refresh
                        sd_image_pigeon.tStopRefresh = tThisFlipGlobal  # on global time
                        sd_image_pigeon.frameNStop = frameN  # exact frame index
                        # update status
                        sd_image_pigeon.status = FINISHED
                        sd_image_pigeon.setAutoDraw(False)
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer], 
                        playbackComponents=[]
                    )
                    # skip the frame we paused on
                    continue
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    simple_demo.forceEnded = routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in simple_demo.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "simple_demo" ---
            for thisComponent in simple_demo.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for simple_demo
            simple_demo.tStop = globalClock.getTime(format='float')
            simple_demo.tStopRefresh = tThisFlipGlobal
            thisExp.addData('simple_demo.stopped', simple_demo.tStop)
            # the Routine "simple_demo" was not non-slip safe, so reset the non-slip timer
            routineTimer.reset()
        # completed 1000.0 repeats of 'simple_demo_trial_loop'
        
    # completed 1.0 repeats of 'simple_demo_block_loop'
    
    
    # set up handler to look after randomisation of conditions etc
    status_instruction_loop = data.TrialHandler2(
        name='status_instruction_loop',
        nReps=1.0, 
        method='sequential', 
        extraInfo=expInfo, 
        originPath=-1, 
        trialList=data.importConditions('Resources/PigeonStatusInstructionsConditions.xlsx'), 
        seed=None, 
    )
    thisExp.addLoop(status_instruction_loop)  # add the loop to the experiment
    thisStatus_instruction_loop = status_instruction_loop.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisStatus_instruction_loop.rgb)
    if thisStatus_instruction_loop != None:
        for paramName in thisStatus_instruction_loop:
            globals()[paramName] = thisStatus_instruction_loop[paramName]
    
    for thisStatus_instruction_loop in status_instruction_loop:
        currentLoop = status_instruction_loop
        thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
        # abbreviate parameter names if possible (e.g. rgb = thisStatus_instruction_loop.rgb)
        if thisStatus_instruction_loop != None:
            for paramName in thisStatus_instruction_loop:
                globals()[paramName] = thisStatus_instruction_loop[paramName]
        
        # --- Prepare to start Routine "status_instruction_setup" ---
        # create an object to store info about Routine status_instruction_setup
        status_instruction_setup = data.Routine(
            name='status_instruction_setup',
            components=[],
        )
        status_instruction_setup.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # Run 'Begin Routine' code from sis_code_setup
        # Set conditions
        status_show = True
        pigeon_show = True
        seeds_show = True
        pigeon_flip = 0.1
        pigeon_true_x = 0
        pigeon_shown_x = pigeon_true_x
        
        box_show = sic_box_show
        box_width = sic_box_width
        box_height = sic_box_height
        box_x = sic_box_x
        box_y = sic_box_y
        message_string = sic_message
        # store start times for status_instruction_setup
        status_instruction_setup.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        status_instruction_setup.tStart = globalClock.getTime(format='float')
        status_instruction_setup.status = STARTED
        thisExp.addData('status_instruction_setup.started', status_instruction_setup.tStart)
        status_instruction_setup.maxDuration = None
        # keep track of which components have finished
        status_instruction_setupComponents = status_instruction_setup.components
        for thisComponent in status_instruction_setup.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "status_instruction_setup" ---
        # if trial has changed, end Routine now
        if isinstance(status_instruction_loop, data.TrialHandler2) and thisStatus_instruction_loop.thisN != status_instruction_loop.thisTrial.thisN:
            continueRoutine = False
        status_instruction_setup.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                status_instruction_setup.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in status_instruction_setup.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "status_instruction_setup" ---
        for thisComponent in status_instruction_setup.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for status_instruction_setup
        status_instruction_setup.tStop = globalClock.getTime(format='float')
        status_instruction_setup.tStopRefresh = tThisFlipGlobal
        thisExp.addData('status_instruction_setup.stopped', status_instruction_setup.tStop)
        # the Routine "status_instruction_setup" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "show_message" ---
        # create an object to store info about Routine show_message
        show_message = data.Routine(
            name='show_message',
            components=[sm_kbd_continue, sm_polygon_box, sm_polygon_abscissa, sm_polygon_ordinate, sm_polygon_petey_midpoint, sm_image_right_seeds, sm_image_left_seeds, sm_image_petey, sm_polygon_coin_bar, sm_polygon_steps_bar, sm_polygon_status_bar_apex, sm_polygon_status_bar_base, sm_text_coins, sm_text_coins_label, sm_text_coins_properties, sm_text_steps, sm_text_steps_label, sm_text_steps_properties, sm_text_message, sm_text_instruction],
        )
        show_message.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # create starting attributes for sm_kbd_continue
        sm_kbd_continue.keys = []
        sm_kbd_continue.rt = []
        _sm_kbd_continue_allKeys = []
        # Run 'Begin Routine' code from sm_code_update
        # Clear all keyboard events
        sm_kbd_continue.clearEvents()
        sm_polygon_box.setPos((box_x,box_y))
        sm_polygon_box.setSize((box_width,box_height))
        sm_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
        sm_polygon_abscissa.setSize((ABSCISSA_WIDTH,LINE_HEIGHT))
        sm_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
        sm_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
        sm_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
        sm_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
        sm_image_petey.setPos((pigeon_shown_x,ORIGIN_Y))
        sm_image_petey.setSize((pigeon_flip, 0.1))
        sm_polygon_coin_bar.setPos((COINS_BAR_X,coins_bar_center))
        sm_polygon_coin_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
        sm_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
        sm_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
        sm_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
        sm_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
        sm_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
        sm_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
        sm_text_coins.setText(coins_count)
        sm_text_coins_properties.setText(coins_string)
        sm_text_steps.setText(steps_max - steps_count)
        sm_text_steps_properties.setText(steps_string)
        sm_text_message.setText(message_string)
        sm_text_instruction.setText(CONTINUE_STRING)
        # store start times for show_message
        show_message.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        show_message.tStart = globalClock.getTime(format='float')
        show_message.status = STARTED
        thisExp.addData('show_message.started', show_message.tStart)
        show_message.maxDuration = None
        # keep track of which components have finished
        show_messageComponents = show_message.components
        for thisComponent in show_message.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "show_message" ---
        # if trial has changed, end Routine now
        if isinstance(status_instruction_loop, data.TrialHandler2) and thisStatus_instruction_loop.thisN != status_instruction_loop.thisTrial.thisN:
            continueRoutine = False
        show_message.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *sm_kbd_continue* updates
            waitOnFlip = False
            
            # if sm_kbd_continue is starting this frame...
            if sm_kbd_continue.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sm_kbd_continue.frameNStart = frameN  # exact frame index
                sm_kbd_continue.tStart = t  # local t and not account for scr refresh
                sm_kbd_continue.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_kbd_continue, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_kbd_continue.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(sm_kbd_continue.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(sm_kbd_continue.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if sm_kbd_continue.status == STARTED and not waitOnFlip:
                theseKeys = sm_kbd_continue.getKeys(keyList=['space', 'up', 'down', 'left', 'right'], ignoreKeys=["escape"], waitRelease=False)
                _sm_kbd_continue_allKeys.extend(theseKeys)
                if len(_sm_kbd_continue_allKeys):
                    sm_kbd_continue.keys = _sm_kbd_continue_allKeys[-1].name  # just the last key pressed
                    sm_kbd_continue.rt = _sm_kbd_continue_allKeys[-1].rt
                    sm_kbd_continue.duration = _sm_kbd_continue_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            
            # *sm_polygon_box* updates
            
            # if sm_polygon_box is starting this frame...
            if sm_polygon_box.status == NOT_STARTED and box_show:
                # keep track of start time/frame for later
                sm_polygon_box.frameNStart = frameN  # exact frame index
                sm_polygon_box.tStart = t  # local t and not account for scr refresh
                sm_polygon_box.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_box, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_box.status = STARTED
                sm_polygon_box.setAutoDraw(True)
            
            # if sm_polygon_box is active this frame...
            if sm_polygon_box.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_abscissa* updates
            
            # if sm_polygon_abscissa is starting this frame...
            if sm_polygon_abscissa.status == NOT_STARTED and pigeon_show:
                # keep track of start time/frame for later
                sm_polygon_abscissa.frameNStart = frameN  # exact frame index
                sm_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                sm_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_abscissa.status = STARTED
                sm_polygon_abscissa.setAutoDraw(True)
            
            # if sm_polygon_abscissa is active this frame...
            if sm_polygon_abscissa.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_ordinate* updates
            
            # if sm_polygon_ordinate is starting this frame...
            if sm_polygon_ordinate.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_ordinate.frameNStart = frameN  # exact frame index
                sm_polygon_ordinate.tStart = t  # local t and not account for scr refresh
                sm_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_ordinate.status = STARTED
                sm_polygon_ordinate.setAutoDraw(True)
            
            # if sm_polygon_ordinate is active this frame...
            if sm_polygon_ordinate.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_petey_midpoint* updates
            
            # if sm_polygon_petey_midpoint is starting this frame...
            if sm_polygon_petey_midpoint.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
                sm_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
                sm_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_petey_midpoint.status = STARTED
                sm_polygon_petey_midpoint.setAutoDraw(True)
            
            # if sm_polygon_petey_midpoint is active this frame...
            if sm_polygon_petey_midpoint.status == STARTED:
                # update params
                pass
            
            # *sm_image_right_seeds* updates
            
            # if sm_image_right_seeds is starting this frame...
            if sm_image_right_seeds.status == NOT_STARTED and seeds_show:
                # keep track of start time/frame for later
                sm_image_right_seeds.frameNStart = frameN  # exact frame index
                sm_image_right_seeds.tStart = t  # local t and not account for scr refresh
                sm_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_right_seeds.status = STARTED
                sm_image_right_seeds.setAutoDraw(True)
            
            # if sm_image_right_seeds is active this frame...
            if sm_image_right_seeds.status == STARTED:
                # update params
                pass
            
            # *sm_image_left_seeds* updates
            
            # if sm_image_left_seeds is starting this frame...
            if sm_image_left_seeds.status == NOT_STARTED and seeds_show:
                # keep track of start time/frame for later
                sm_image_left_seeds.frameNStart = frameN  # exact frame index
                sm_image_left_seeds.tStart = t  # local t and not account for scr refresh
                sm_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_left_seeds.status = STARTED
                sm_image_left_seeds.setAutoDraw(True)
            
            # if sm_image_left_seeds is active this frame...
            if sm_image_left_seeds.status == STARTED:
                # update params
                pass
            
            # *sm_image_petey* updates
            
            # if sm_image_petey is starting this frame...
            if sm_image_petey.status == NOT_STARTED and pigeon_show:
                # keep track of start time/frame for later
                sm_image_petey.frameNStart = frameN  # exact frame index
                sm_image_petey.tStart = t  # local t and not account for scr refresh
                sm_image_petey.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_petey, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_petey.status = STARTED
                sm_image_petey.setAutoDraw(True)
            
            # if sm_image_petey is active this frame...
            if sm_image_petey.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_coin_bar* updates
            
            # if sm_polygon_coin_bar is starting this frame...
            if sm_polygon_coin_bar.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_coin_bar.frameNStart = frameN  # exact frame index
                sm_polygon_coin_bar.tStart = t  # local t and not account for scr refresh
                sm_polygon_coin_bar.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_coin_bar, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_coin_bar.status = STARTED
                sm_polygon_coin_bar.setAutoDraw(True)
            
            # if sm_polygon_coin_bar is active this frame...
            if sm_polygon_coin_bar.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_steps_bar* updates
            
            # if sm_polygon_steps_bar is starting this frame...
            if sm_polygon_steps_bar.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_steps_bar.frameNStart = frameN  # exact frame index
                sm_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
                sm_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_steps_bar.status = STARTED
                sm_polygon_steps_bar.setAutoDraw(True)
            
            # if sm_polygon_steps_bar is active this frame...
            if sm_polygon_steps_bar.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_status_bar_apex* updates
            
            # if sm_polygon_status_bar_apex is starting this frame...
            if sm_polygon_status_bar_apex.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
                sm_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
                sm_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_status_bar_apex.status = STARTED
                sm_polygon_status_bar_apex.setAutoDraw(True)
            
            # if sm_polygon_status_bar_apex is active this frame...
            if sm_polygon_status_bar_apex.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_status_bar_base* updates
            
            # if sm_polygon_status_bar_base is starting this frame...
            if sm_polygon_status_bar_base.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_status_bar_base.frameNStart = frameN  # exact frame index
                sm_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
                sm_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_status_bar_base.status = STARTED
                sm_polygon_status_bar_base.setAutoDraw(True)
            
            # if sm_polygon_status_bar_base is active this frame...
            if sm_polygon_status_bar_base.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins* updates
            
            # if sm_text_coins is starting this frame...
            if sm_text_coins.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins.frameNStart = frameN  # exact frame index
                sm_text_coins.tStart = t  # local t and not account for scr refresh
                sm_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins.status = STARTED
                sm_text_coins.setAutoDraw(True)
            
            # if sm_text_coins is active this frame...
            if sm_text_coins.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins_label* updates
            
            # if sm_text_coins_label is starting this frame...
            if sm_text_coins_label.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins_label.frameNStart = frameN  # exact frame index
                sm_text_coins_label.tStart = t  # local t and not account for scr refresh
                sm_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins_label, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins_label.status = STARTED
                sm_text_coins_label.setAutoDraw(True)
            
            # if sm_text_coins_label is active this frame...
            if sm_text_coins_label.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins_properties* updates
            
            # if sm_text_coins_properties is starting this frame...
            if sm_text_coins_properties.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins_properties.frameNStart = frameN  # exact frame index
                sm_text_coins_properties.tStart = t  # local t and not account for scr refresh
                sm_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins_properties.status = STARTED
                sm_text_coins_properties.setAutoDraw(True)
            
            # if sm_text_coins_properties is active this frame...
            if sm_text_coins_properties.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps* updates
            
            # if sm_text_steps is starting this frame...
            if sm_text_steps.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps.frameNStart = frameN  # exact frame index
                sm_text_steps.tStart = t  # local t and not account for scr refresh
                sm_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps.status = STARTED
                sm_text_steps.setAutoDraw(True)
            
            # if sm_text_steps is active this frame...
            if sm_text_steps.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps_label* updates
            
            # if sm_text_steps_label is starting this frame...
            if sm_text_steps_label.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps_label.frameNStart = frameN  # exact frame index
                sm_text_steps_label.tStart = t  # local t and not account for scr refresh
                sm_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps_label, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps_label.status = STARTED
                sm_text_steps_label.setAutoDraw(True)
            
            # if sm_text_steps_label is active this frame...
            if sm_text_steps_label.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps_properties* updates
            
            # if sm_text_steps_properties is starting this frame...
            if sm_text_steps_properties.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps_properties.frameNStart = frameN  # exact frame index
                sm_text_steps_properties.tStart = t  # local t and not account for scr refresh
                sm_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps_properties.status = STARTED
                sm_text_steps_properties.setAutoDraw(True)
            
            # if sm_text_steps_properties is active this frame...
            if sm_text_steps_properties.status == STARTED:
                # update params
                pass
            
            # *sm_text_message* updates
            
            # if sm_text_message is starting this frame...
            if sm_text_message.status == NOT_STARTED and message_show:
                # keep track of start time/frame for later
                sm_text_message.frameNStart = frameN  # exact frame index
                sm_text_message.tStart = t  # local t and not account for scr refresh
                sm_text_message.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_message, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_message.status = STARTED
                sm_text_message.setAutoDraw(True)
            
            # if sm_text_message is active this frame...
            if sm_text_message.status == STARTED:
                # update params
                pass
            
            # *sm_text_instruction* updates
            
            # if sm_text_instruction is starting this frame...
            if sm_text_instruction.status == NOT_STARTED and instruction_show:
                # keep track of start time/frame for later
                sm_text_instruction.frameNStart = frameN  # exact frame index
                sm_text_instruction.tStart = t  # local t and not account for scr refresh
                sm_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_instruction, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_instruction.status = STARTED
                sm_text_instruction.setAutoDraw(True)
            
            # if sm_text_instruction is active this frame...
            if sm_text_instruction.status == STARTED:
                # update params
                pass
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                show_message.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in show_message.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "show_message" ---
        for thisComponent in show_message.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for show_message
        show_message.tStop = globalClock.getTime(format='float')
        show_message.tStopRefresh = tThisFlipGlobal
        thisExp.addData('show_message.stopped', show_message.tStop)
        # Run 'End Routine' code from sm_code_update
        # Ick, but whatever. Makes sure that we don't ever
        # drop unexpectedly out of a move loop without 
        # resetting
        trial_start = True
        # the Routine "show_message" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
    # completed 1.0 repeats of 'status_instruction_loop'
    
    
    # set up handler to look after randomisation of conditions etc
    experiment_loop = data.TrialHandler2(
        name='experiment_loop',
        nReps=1.0, 
        method='sequential', 
        extraInfo=expInfo, 
        originPath=-1, 
        trialList=data.importConditions('Resources/PigeonTaskConditions.xlsx'), 
        seed=None, 
    )
    thisExp.addLoop(experiment_loop)  # add the loop to the experiment
    thisExperiment_loop = experiment_loop.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisExperiment_loop.rgb)
    if thisExperiment_loop != None:
        for paramName in thisExperiment_loop:
            globals()[paramName] = thisExperiment_loop[paramName]
    if thisSession is not None:
        # if running in a Session with a Liaison client, send data up to now
        thisSession.sendExperimentData()
    
    for thisExperiment_loop in experiment_loop:
        currentLoop = experiment_loop
        thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        # abbreviate parameter names if possible (e.g. rgb = thisExperiment_loop.rgb)
        if thisExperiment_loop != None:
            for paramName in thisExperiment_loop:
                globals()[paramName] = thisExperiment_loop[paramName]
        
        # --- Prepare to start Routine "start_task" ---
        # create an object to store info about Routine start_task
        start_task = data.Routine(
            name='start_task',
            components=[],
        )
        start_task.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # Run 'Begin Routine' code from st_code_setup
        # Setup task
        trial_start = True
        trial_number = 1
        steps_count = 0
        coins_count = 0
        last_coins_count = coins_count
        last_steps_count = steps_count
        
        # Set globals from conditions file variables
        task_type = tc_task_type
        block_name = tc_block_name
        steps_max = tc_steps_max
        steps_taken_to_start_trial = tc_steps_taken_to_start_trial
        steps_lost_per_error = tc_steps_lost_per_error
        coins_max = tc_coins_max
        coins_paid_to_start_trial = tc_coins_paid_to_start_trial
        coins_lost_per_error = tc_coins_lost_per_error
        coins_gained_per_correct = tc_coins_gained_per_correct
        predefined_bound_min_bonus = tc_predefined_bound_min_bonus 
        predefined_bound_max_bonus = tc_predefined_bound_max_bonus
        snr_set = tc_snr_set
        
        # Get SNR conditions file name
        if snr_set < 10:
            snr_conditions_file = 'Resources/PigeonSNRSet0' + str(snr_set) + 'Conditions.xlsx'
        else:
            snr_conditions_file = 'Resources/PigeonSNRSet' + str(snr_set) + 'Conditions.xlsx'    
        
        # Get Instructions Conditions file name
        instructions_conditions_file = 'Resources/Pigeon_' + task_type + '_InstructionsConditions.xlsx'
        
        # Control flow
        if task_type == 'online':
            online_trials_count = 10000
            predefined_trials_count = 0
        else:
            online_trials_count = 0
            predefined_trials_count = 10000
            predefined_bound_x = PREDEFINED_BOUND_X_DEFAULT
            
        # New block
        block_number += 1
        message_string = 'Block ' + str(block_number) + ' of ' + str(len(experiment_loop.trialList)) + '\n\n' + 'Type = ' + task_type
        
        # Initialize pigeon variables
        pigeon_flip = 0.1
        pigeon_true_x = 0
        pigeon_shown_x = pigeon_true_x
        
        # Update status graphics/text
        update_status(True, True, True)
        # store start times for start_task
        start_task.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        start_task.tStart = globalClock.getTime(format='float')
        start_task.status = STARTED
        thisExp.addData('start_task.started', start_task.tStart)
        start_task.maxDuration = None
        # keep track of which components have finished
        start_taskComponents = start_task.components
        for thisComponent in start_task.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "start_task" ---
        # if trial has changed, end Routine now
        if isinstance(experiment_loop, data.TrialHandler2) and thisExperiment_loop.thisN != experiment_loop.thisTrial.thisN:
            continueRoutine = False
        start_task.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                start_task.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in start_task.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "start_task" ---
        for thisComponent in start_task.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for start_task
        start_task.tStop = globalClock.getTime(format='float')
        start_task.tStopRefresh = tThisFlipGlobal
        thisExp.addData('start_task.stopped', start_task.tStop)
        # the Routine "start_task" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "show_message" ---
        # create an object to store info about Routine show_message
        show_message = data.Routine(
            name='show_message',
            components=[sm_kbd_continue, sm_polygon_box, sm_polygon_abscissa, sm_polygon_ordinate, sm_polygon_petey_midpoint, sm_image_right_seeds, sm_image_left_seeds, sm_image_petey, sm_polygon_coin_bar, sm_polygon_steps_bar, sm_polygon_status_bar_apex, sm_polygon_status_bar_base, sm_text_coins, sm_text_coins_label, sm_text_coins_properties, sm_text_steps, sm_text_steps_label, sm_text_steps_properties, sm_text_message, sm_text_instruction],
        )
        show_message.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # create starting attributes for sm_kbd_continue
        sm_kbd_continue.keys = []
        sm_kbd_continue.rt = []
        _sm_kbd_continue_allKeys = []
        # Run 'Begin Routine' code from sm_code_update
        # Clear all keyboard events
        sm_kbd_continue.clearEvents()
        sm_polygon_box.setPos((box_x,box_y))
        sm_polygon_box.setSize((box_width,box_height))
        sm_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
        sm_polygon_abscissa.setSize((ABSCISSA_WIDTH,LINE_HEIGHT))
        sm_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
        sm_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
        sm_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
        sm_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
        sm_image_petey.setPos((pigeon_shown_x,ORIGIN_Y))
        sm_image_petey.setSize((pigeon_flip, 0.1))
        sm_polygon_coin_bar.setPos((COINS_BAR_X,coins_bar_center))
        sm_polygon_coin_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
        sm_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
        sm_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
        sm_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
        sm_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
        sm_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
        sm_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
        sm_text_coins.setText(coins_count)
        sm_text_coins_properties.setText(coins_string)
        sm_text_steps.setText(steps_max - steps_count)
        sm_text_steps_properties.setText(steps_string)
        sm_text_message.setText(message_string)
        sm_text_instruction.setText(CONTINUE_STRING)
        # store start times for show_message
        show_message.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        show_message.tStart = globalClock.getTime(format='float')
        show_message.status = STARTED
        thisExp.addData('show_message.started', show_message.tStart)
        show_message.maxDuration = None
        # keep track of which components have finished
        show_messageComponents = show_message.components
        for thisComponent in show_message.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "show_message" ---
        # if trial has changed, end Routine now
        if isinstance(experiment_loop, data.TrialHandler2) and thisExperiment_loop.thisN != experiment_loop.thisTrial.thisN:
            continueRoutine = False
        show_message.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *sm_kbd_continue* updates
            waitOnFlip = False
            
            # if sm_kbd_continue is starting this frame...
            if sm_kbd_continue.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sm_kbd_continue.frameNStart = frameN  # exact frame index
                sm_kbd_continue.tStart = t  # local t and not account for scr refresh
                sm_kbd_continue.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_kbd_continue, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_kbd_continue.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(sm_kbd_continue.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(sm_kbd_continue.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if sm_kbd_continue.status == STARTED and not waitOnFlip:
                theseKeys = sm_kbd_continue.getKeys(keyList=['space', 'up', 'down', 'left', 'right'], ignoreKeys=["escape"], waitRelease=False)
                _sm_kbd_continue_allKeys.extend(theseKeys)
                if len(_sm_kbd_continue_allKeys):
                    sm_kbd_continue.keys = _sm_kbd_continue_allKeys[-1].name  # just the last key pressed
                    sm_kbd_continue.rt = _sm_kbd_continue_allKeys[-1].rt
                    sm_kbd_continue.duration = _sm_kbd_continue_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            
            # *sm_polygon_box* updates
            
            # if sm_polygon_box is starting this frame...
            if sm_polygon_box.status == NOT_STARTED and box_show:
                # keep track of start time/frame for later
                sm_polygon_box.frameNStart = frameN  # exact frame index
                sm_polygon_box.tStart = t  # local t and not account for scr refresh
                sm_polygon_box.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_box, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_box.status = STARTED
                sm_polygon_box.setAutoDraw(True)
            
            # if sm_polygon_box is active this frame...
            if sm_polygon_box.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_abscissa* updates
            
            # if sm_polygon_abscissa is starting this frame...
            if sm_polygon_abscissa.status == NOT_STARTED and pigeon_show:
                # keep track of start time/frame for later
                sm_polygon_abscissa.frameNStart = frameN  # exact frame index
                sm_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                sm_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_abscissa.status = STARTED
                sm_polygon_abscissa.setAutoDraw(True)
            
            # if sm_polygon_abscissa is active this frame...
            if sm_polygon_abscissa.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_ordinate* updates
            
            # if sm_polygon_ordinate is starting this frame...
            if sm_polygon_ordinate.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_ordinate.frameNStart = frameN  # exact frame index
                sm_polygon_ordinate.tStart = t  # local t and not account for scr refresh
                sm_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_ordinate.status = STARTED
                sm_polygon_ordinate.setAutoDraw(True)
            
            # if sm_polygon_ordinate is active this frame...
            if sm_polygon_ordinate.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_petey_midpoint* updates
            
            # if sm_polygon_petey_midpoint is starting this frame...
            if sm_polygon_petey_midpoint.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
                sm_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
                sm_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_petey_midpoint.status = STARTED
                sm_polygon_petey_midpoint.setAutoDraw(True)
            
            # if sm_polygon_petey_midpoint is active this frame...
            if sm_polygon_petey_midpoint.status == STARTED:
                # update params
                pass
            
            # *sm_image_right_seeds* updates
            
            # if sm_image_right_seeds is starting this frame...
            if sm_image_right_seeds.status == NOT_STARTED and seeds_show:
                # keep track of start time/frame for later
                sm_image_right_seeds.frameNStart = frameN  # exact frame index
                sm_image_right_seeds.tStart = t  # local t and not account for scr refresh
                sm_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_right_seeds.status = STARTED
                sm_image_right_seeds.setAutoDraw(True)
            
            # if sm_image_right_seeds is active this frame...
            if sm_image_right_seeds.status == STARTED:
                # update params
                pass
            
            # *sm_image_left_seeds* updates
            
            # if sm_image_left_seeds is starting this frame...
            if sm_image_left_seeds.status == NOT_STARTED and seeds_show:
                # keep track of start time/frame for later
                sm_image_left_seeds.frameNStart = frameN  # exact frame index
                sm_image_left_seeds.tStart = t  # local t and not account for scr refresh
                sm_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_left_seeds.status = STARTED
                sm_image_left_seeds.setAutoDraw(True)
            
            # if sm_image_left_seeds is active this frame...
            if sm_image_left_seeds.status == STARTED:
                # update params
                pass
            
            # *sm_image_petey* updates
            
            # if sm_image_petey is starting this frame...
            if sm_image_petey.status == NOT_STARTED and pigeon_show:
                # keep track of start time/frame for later
                sm_image_petey.frameNStart = frameN  # exact frame index
                sm_image_petey.tStart = t  # local t and not account for scr refresh
                sm_image_petey.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_image_petey, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_image_petey.status = STARTED
                sm_image_petey.setAutoDraw(True)
            
            # if sm_image_petey is active this frame...
            if sm_image_petey.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_coin_bar* updates
            
            # if sm_polygon_coin_bar is starting this frame...
            if sm_polygon_coin_bar.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_coin_bar.frameNStart = frameN  # exact frame index
                sm_polygon_coin_bar.tStart = t  # local t and not account for scr refresh
                sm_polygon_coin_bar.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_coin_bar, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_coin_bar.status = STARTED
                sm_polygon_coin_bar.setAutoDraw(True)
            
            # if sm_polygon_coin_bar is active this frame...
            if sm_polygon_coin_bar.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_steps_bar* updates
            
            # if sm_polygon_steps_bar is starting this frame...
            if sm_polygon_steps_bar.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_steps_bar.frameNStart = frameN  # exact frame index
                sm_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
                sm_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_steps_bar.status = STARTED
                sm_polygon_steps_bar.setAutoDraw(True)
            
            # if sm_polygon_steps_bar is active this frame...
            if sm_polygon_steps_bar.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_status_bar_apex* updates
            
            # if sm_polygon_status_bar_apex is starting this frame...
            if sm_polygon_status_bar_apex.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
                sm_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
                sm_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_status_bar_apex.status = STARTED
                sm_polygon_status_bar_apex.setAutoDraw(True)
            
            # if sm_polygon_status_bar_apex is active this frame...
            if sm_polygon_status_bar_apex.status == STARTED:
                # update params
                pass
            
            # *sm_polygon_status_bar_base* updates
            
            # if sm_polygon_status_bar_base is starting this frame...
            if sm_polygon_status_bar_base.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_polygon_status_bar_base.frameNStart = frameN  # exact frame index
                sm_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
                sm_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_polygon_status_bar_base.status = STARTED
                sm_polygon_status_bar_base.setAutoDraw(True)
            
            # if sm_polygon_status_bar_base is active this frame...
            if sm_polygon_status_bar_base.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins* updates
            
            # if sm_text_coins is starting this frame...
            if sm_text_coins.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins.frameNStart = frameN  # exact frame index
                sm_text_coins.tStart = t  # local t and not account for scr refresh
                sm_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins.status = STARTED
                sm_text_coins.setAutoDraw(True)
            
            # if sm_text_coins is active this frame...
            if sm_text_coins.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins_label* updates
            
            # if sm_text_coins_label is starting this frame...
            if sm_text_coins_label.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins_label.frameNStart = frameN  # exact frame index
                sm_text_coins_label.tStart = t  # local t and not account for scr refresh
                sm_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins_label, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins_label.status = STARTED
                sm_text_coins_label.setAutoDraw(True)
            
            # if sm_text_coins_label is active this frame...
            if sm_text_coins_label.status == STARTED:
                # update params
                pass
            
            # *sm_text_coins_properties* updates
            
            # if sm_text_coins_properties is starting this frame...
            if sm_text_coins_properties.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_coins_properties.frameNStart = frameN  # exact frame index
                sm_text_coins_properties.tStart = t  # local t and not account for scr refresh
                sm_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_coins_properties.status = STARTED
                sm_text_coins_properties.setAutoDraw(True)
            
            # if sm_text_coins_properties is active this frame...
            if sm_text_coins_properties.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps* updates
            
            # if sm_text_steps is starting this frame...
            if sm_text_steps.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps.frameNStart = frameN  # exact frame index
                sm_text_steps.tStart = t  # local t and not account for scr refresh
                sm_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps.status = STARTED
                sm_text_steps.setAutoDraw(True)
            
            # if sm_text_steps is active this frame...
            if sm_text_steps.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps_label* updates
            
            # if sm_text_steps_label is starting this frame...
            if sm_text_steps_label.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps_label.frameNStart = frameN  # exact frame index
                sm_text_steps_label.tStart = t  # local t and not account for scr refresh
                sm_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps_label, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps_label.status = STARTED
                sm_text_steps_label.setAutoDraw(True)
            
            # if sm_text_steps_label is active this frame...
            if sm_text_steps_label.status == STARTED:
                # update params
                pass
            
            # *sm_text_steps_properties* updates
            
            # if sm_text_steps_properties is starting this frame...
            if sm_text_steps_properties.status == NOT_STARTED and status_show:
                # keep track of start time/frame for later
                sm_text_steps_properties.frameNStart = frameN  # exact frame index
                sm_text_steps_properties.tStart = t  # local t and not account for scr refresh
                sm_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_steps_properties.status = STARTED
                sm_text_steps_properties.setAutoDraw(True)
            
            # if sm_text_steps_properties is active this frame...
            if sm_text_steps_properties.status == STARTED:
                # update params
                pass
            
            # *sm_text_message* updates
            
            # if sm_text_message is starting this frame...
            if sm_text_message.status == NOT_STARTED and message_show:
                # keep track of start time/frame for later
                sm_text_message.frameNStart = frameN  # exact frame index
                sm_text_message.tStart = t  # local t and not account for scr refresh
                sm_text_message.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_message, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_message.status = STARTED
                sm_text_message.setAutoDraw(True)
            
            # if sm_text_message is active this frame...
            if sm_text_message.status == STARTED:
                # update params
                pass
            
            # *sm_text_instruction* updates
            
            # if sm_text_instruction is starting this frame...
            if sm_text_instruction.status == NOT_STARTED and instruction_show:
                # keep track of start time/frame for later
                sm_text_instruction.frameNStart = frameN  # exact frame index
                sm_text_instruction.tStart = t  # local t and not account for scr refresh
                sm_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sm_text_instruction, 'tStartRefresh')  # time at next scr refresh
                # update status
                sm_text_instruction.status = STARTED
                sm_text_instruction.setAutoDraw(True)
            
            # if sm_text_instruction is active this frame...
            if sm_text_instruction.status == STARTED:
                # update params
                pass
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer], 
                    playbackComponents=[]
                )
                # skip the frame we paused on
                continue
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                show_message.forceEnded = routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in show_message.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "show_message" ---
        for thisComponent in show_message.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for show_message
        show_message.tStop = globalClock.getTime(format='float')
        show_message.tStopRefresh = tThisFlipGlobal
        thisExp.addData('show_message.stopped', show_message.tStop)
        # Run 'End Routine' code from sm_code_update
        # Ick, but whatever. Makes sure that we don't ever
        # drop unexpectedly out of a move loop without 
        # resetting
        trial_start = True
        # the Routine "show_message" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # set up handler to look after randomisation of conditions etc
        instructions_loop = data.TrialHandler2(
            name='instructions_loop',
            nReps=1.0, 
            method='sequential', 
            extraInfo=expInfo, 
            originPath=-1, 
            trialList=data.importConditions(instructions_conditions_file), 
            seed=None, 
        )
        thisExp.addLoop(instructions_loop)  # add the loop to the experiment
        thisInstructions_loop = instructions_loop.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb = thisInstructions_loop.rgb)
        if thisInstructions_loop != None:
            for paramName in thisInstructions_loop:
                globals()[paramName] = thisInstructions_loop[paramName]
        
        for thisInstructions_loop in instructions_loop:
            currentLoop = instructions_loop
            thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
            # abbreviate parameter names if possible (e.g. rgb = thisInstructions_loop.rgb)
            if thisInstructions_loop != None:
                for paramName in thisInstructions_loop:
                    globals()[paramName] = thisInstructions_loop[paramName]
            
            # --- Prepare to start Routine "setup_instructions_message" ---
            # create an object to store info about Routine setup_instructions_message
            setup_instructions_message = data.Routine(
                name='setup_instructions_message',
                components=[],
            )
            setup_instructions_message.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            # Run 'Begin Routine' code from code
            message_string = ic_message_string
            # store start times for setup_instructions_message
            setup_instructions_message.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            setup_instructions_message.tStart = globalClock.getTime(format='float')
            setup_instructions_message.status = STARTED
            thisExp.addData('setup_instructions_message.started', setup_instructions_message.tStart)
            setup_instructions_message.maxDuration = None
            # keep track of which components have finished
            setup_instructions_messageComponents = setup_instructions_message.components
            for thisComponent in setup_instructions_message.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "setup_instructions_message" ---
            # if trial has changed, end Routine now
            if isinstance(instructions_loop, data.TrialHandler2) and thisInstructions_loop.thisN != instructions_loop.thisTrial.thisN:
                continueRoutine = False
            setup_instructions_message.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer], 
                        playbackComponents=[]
                    )
                    # skip the frame we paused on
                    continue
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    setup_instructions_message.forceEnded = routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in setup_instructions_message.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "setup_instructions_message" ---
            for thisComponent in setup_instructions_message.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for setup_instructions_message
            setup_instructions_message.tStop = globalClock.getTime(format='float')
            setup_instructions_message.tStopRefresh = tThisFlipGlobal
            thisExp.addData('setup_instructions_message.stopped', setup_instructions_message.tStop)
            # the Routine "setup_instructions_message" was not non-slip safe, so reset the non-slip timer
            routineTimer.reset()
            
            # --- Prepare to start Routine "show_message" ---
            # create an object to store info about Routine show_message
            show_message = data.Routine(
                name='show_message',
                components=[sm_kbd_continue, sm_polygon_box, sm_polygon_abscissa, sm_polygon_ordinate, sm_polygon_petey_midpoint, sm_image_right_seeds, sm_image_left_seeds, sm_image_petey, sm_polygon_coin_bar, sm_polygon_steps_bar, sm_polygon_status_bar_apex, sm_polygon_status_bar_base, sm_text_coins, sm_text_coins_label, sm_text_coins_properties, sm_text_steps, sm_text_steps_label, sm_text_steps_properties, sm_text_message, sm_text_instruction],
            )
            show_message.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            # create starting attributes for sm_kbd_continue
            sm_kbd_continue.keys = []
            sm_kbd_continue.rt = []
            _sm_kbd_continue_allKeys = []
            # Run 'Begin Routine' code from sm_code_update
            # Clear all keyboard events
            sm_kbd_continue.clearEvents()
            sm_polygon_box.setPos((box_x,box_y))
            sm_polygon_box.setSize((box_width,box_height))
            sm_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
            sm_polygon_abscissa.setSize((ABSCISSA_WIDTH,LINE_HEIGHT))
            sm_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
            sm_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
            sm_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
            sm_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
            sm_image_petey.setPos((pigeon_shown_x,ORIGIN_Y))
            sm_image_petey.setSize((pigeon_flip, 0.1))
            sm_polygon_coin_bar.setPos((COINS_BAR_X,coins_bar_center))
            sm_polygon_coin_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
            sm_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
            sm_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
            sm_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
            sm_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
            sm_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
            sm_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
            sm_text_coins.setText(coins_count)
            sm_text_coins_properties.setText(coins_string)
            sm_text_steps.setText(steps_max - steps_count)
            sm_text_steps_properties.setText(steps_string)
            sm_text_message.setText(message_string)
            sm_text_instruction.setText(CONTINUE_STRING)
            # store start times for show_message
            show_message.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            show_message.tStart = globalClock.getTime(format='float')
            show_message.status = STARTED
            thisExp.addData('show_message.started', show_message.tStart)
            show_message.maxDuration = None
            # keep track of which components have finished
            show_messageComponents = show_message.components
            for thisComponent in show_message.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "show_message" ---
            # if trial has changed, end Routine now
            if isinstance(instructions_loop, data.TrialHandler2) and thisInstructions_loop.thisN != instructions_loop.thisTrial.thisN:
                continueRoutine = False
            show_message.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                # *sm_kbd_continue* updates
                waitOnFlip = False
                
                # if sm_kbd_continue is starting this frame...
                if sm_kbd_continue.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    sm_kbd_continue.frameNStart = frameN  # exact frame index
                    sm_kbd_continue.tStart = t  # local t and not account for scr refresh
                    sm_kbd_continue.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_kbd_continue, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_kbd_continue.status = STARTED
                    # keyboard checking is just starting
                    waitOnFlip = True
                    win.callOnFlip(sm_kbd_continue.clock.reset)  # t=0 on next screen flip
                    win.callOnFlip(sm_kbd_continue.clearEvents, eventType='keyboard')  # clear events on next screen flip
                if sm_kbd_continue.status == STARTED and not waitOnFlip:
                    theseKeys = sm_kbd_continue.getKeys(keyList=['space', 'up', 'down', 'left', 'right'], ignoreKeys=["escape"], waitRelease=False)
                    _sm_kbd_continue_allKeys.extend(theseKeys)
                    if len(_sm_kbd_continue_allKeys):
                        sm_kbd_continue.keys = _sm_kbd_continue_allKeys[-1].name  # just the last key pressed
                        sm_kbd_continue.rt = _sm_kbd_continue_allKeys[-1].rt
                        sm_kbd_continue.duration = _sm_kbd_continue_allKeys[-1].duration
                        # a response ends the routine
                        continueRoutine = False
                
                # *sm_polygon_box* updates
                
                # if sm_polygon_box is starting this frame...
                if sm_polygon_box.status == NOT_STARTED and box_show:
                    # keep track of start time/frame for later
                    sm_polygon_box.frameNStart = frameN  # exact frame index
                    sm_polygon_box.tStart = t  # local t and not account for scr refresh
                    sm_polygon_box.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_box, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_box.status = STARTED
                    sm_polygon_box.setAutoDraw(True)
                
                # if sm_polygon_box is active this frame...
                if sm_polygon_box.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_abscissa* updates
                
                # if sm_polygon_abscissa is starting this frame...
                if sm_polygon_abscissa.status == NOT_STARTED and pigeon_show:
                    # keep track of start time/frame for later
                    sm_polygon_abscissa.frameNStart = frameN  # exact frame index
                    sm_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                    sm_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_abscissa.status = STARTED
                    sm_polygon_abscissa.setAutoDraw(True)
                
                # if sm_polygon_abscissa is active this frame...
                if sm_polygon_abscissa.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_ordinate* updates
                
                # if sm_polygon_ordinate is starting this frame...
                if sm_polygon_ordinate.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_ordinate.frameNStart = frameN  # exact frame index
                    sm_polygon_ordinate.tStart = t  # local t and not account for scr refresh
                    sm_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_ordinate.status = STARTED
                    sm_polygon_ordinate.setAutoDraw(True)
                
                # if sm_polygon_ordinate is active this frame...
                if sm_polygon_ordinate.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_petey_midpoint* updates
                
                # if sm_polygon_petey_midpoint is starting this frame...
                if sm_polygon_petey_midpoint.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
                    sm_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
                    sm_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_petey_midpoint.status = STARTED
                    sm_polygon_petey_midpoint.setAutoDraw(True)
                
                # if sm_polygon_petey_midpoint is active this frame...
                if sm_polygon_petey_midpoint.status == STARTED:
                    # update params
                    pass
                
                # *sm_image_right_seeds* updates
                
                # if sm_image_right_seeds is starting this frame...
                if sm_image_right_seeds.status == NOT_STARTED and seeds_show:
                    # keep track of start time/frame for later
                    sm_image_right_seeds.frameNStart = frameN  # exact frame index
                    sm_image_right_seeds.tStart = t  # local t and not account for scr refresh
                    sm_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_image_right_seeds.status = STARTED
                    sm_image_right_seeds.setAutoDraw(True)
                
                # if sm_image_right_seeds is active this frame...
                if sm_image_right_seeds.status == STARTED:
                    # update params
                    pass
                
                # *sm_image_left_seeds* updates
                
                # if sm_image_left_seeds is starting this frame...
                if sm_image_left_seeds.status == NOT_STARTED and seeds_show:
                    # keep track of start time/frame for later
                    sm_image_left_seeds.frameNStart = frameN  # exact frame index
                    sm_image_left_seeds.tStart = t  # local t and not account for scr refresh
                    sm_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_image_left_seeds.status = STARTED
                    sm_image_left_seeds.setAutoDraw(True)
                
                # if sm_image_left_seeds is active this frame...
                if sm_image_left_seeds.status == STARTED:
                    # update params
                    pass
                
                # *sm_image_petey* updates
                
                # if sm_image_petey is starting this frame...
                if sm_image_petey.status == NOT_STARTED and pigeon_show:
                    # keep track of start time/frame for later
                    sm_image_petey.frameNStart = frameN  # exact frame index
                    sm_image_petey.tStart = t  # local t and not account for scr refresh
                    sm_image_petey.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_image_petey, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_image_petey.status = STARTED
                    sm_image_petey.setAutoDraw(True)
                
                # if sm_image_petey is active this frame...
                if sm_image_petey.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_coin_bar* updates
                
                # if sm_polygon_coin_bar is starting this frame...
                if sm_polygon_coin_bar.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_coin_bar.frameNStart = frameN  # exact frame index
                    sm_polygon_coin_bar.tStart = t  # local t and not account for scr refresh
                    sm_polygon_coin_bar.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_coin_bar, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_coin_bar.status = STARTED
                    sm_polygon_coin_bar.setAutoDraw(True)
                
                # if sm_polygon_coin_bar is active this frame...
                if sm_polygon_coin_bar.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_steps_bar* updates
                
                # if sm_polygon_steps_bar is starting this frame...
                if sm_polygon_steps_bar.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_steps_bar.frameNStart = frameN  # exact frame index
                    sm_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
                    sm_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_steps_bar.status = STARTED
                    sm_polygon_steps_bar.setAutoDraw(True)
                
                # if sm_polygon_steps_bar is active this frame...
                if sm_polygon_steps_bar.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_status_bar_apex* updates
                
                # if sm_polygon_status_bar_apex is starting this frame...
                if sm_polygon_status_bar_apex.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
                    sm_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
                    sm_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_status_bar_apex.status = STARTED
                    sm_polygon_status_bar_apex.setAutoDraw(True)
                
                # if sm_polygon_status_bar_apex is active this frame...
                if sm_polygon_status_bar_apex.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_status_bar_base* updates
                
                # if sm_polygon_status_bar_base is starting this frame...
                if sm_polygon_status_bar_base.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_status_bar_base.frameNStart = frameN  # exact frame index
                    sm_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
                    sm_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_status_bar_base.status = STARTED
                    sm_polygon_status_bar_base.setAutoDraw(True)
                
                # if sm_polygon_status_bar_base is active this frame...
                if sm_polygon_status_bar_base.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_coins* updates
                
                # if sm_text_coins is starting this frame...
                if sm_text_coins.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_coins.frameNStart = frameN  # exact frame index
                    sm_text_coins.tStart = t  # local t and not account for scr refresh
                    sm_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_coins, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_coins.status = STARTED
                    sm_text_coins.setAutoDraw(True)
                
                # if sm_text_coins is active this frame...
                if sm_text_coins.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_coins_label* updates
                
                # if sm_text_coins_label is starting this frame...
                if sm_text_coins_label.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_coins_label.frameNStart = frameN  # exact frame index
                    sm_text_coins_label.tStart = t  # local t and not account for scr refresh
                    sm_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_coins_label, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_coins_label.status = STARTED
                    sm_text_coins_label.setAutoDraw(True)
                
                # if sm_text_coins_label is active this frame...
                if sm_text_coins_label.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_coins_properties* updates
                
                # if sm_text_coins_properties is starting this frame...
                if sm_text_coins_properties.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_coins_properties.frameNStart = frameN  # exact frame index
                    sm_text_coins_properties.tStart = t  # local t and not account for scr refresh
                    sm_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_coins_properties.status = STARTED
                    sm_text_coins_properties.setAutoDraw(True)
                
                # if sm_text_coins_properties is active this frame...
                if sm_text_coins_properties.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_steps* updates
                
                # if sm_text_steps is starting this frame...
                if sm_text_steps.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_steps.frameNStart = frameN  # exact frame index
                    sm_text_steps.tStart = t  # local t and not account for scr refresh
                    sm_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_steps, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_steps.status = STARTED
                    sm_text_steps.setAutoDraw(True)
                
                # if sm_text_steps is active this frame...
                if sm_text_steps.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_steps_label* updates
                
                # if sm_text_steps_label is starting this frame...
                if sm_text_steps_label.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_steps_label.frameNStart = frameN  # exact frame index
                    sm_text_steps_label.tStart = t  # local t and not account for scr refresh
                    sm_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_steps_label, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_steps_label.status = STARTED
                    sm_text_steps_label.setAutoDraw(True)
                
                # if sm_text_steps_label is active this frame...
                if sm_text_steps_label.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_steps_properties* updates
                
                # if sm_text_steps_properties is starting this frame...
                if sm_text_steps_properties.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_steps_properties.frameNStart = frameN  # exact frame index
                    sm_text_steps_properties.tStart = t  # local t and not account for scr refresh
                    sm_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_steps_properties.status = STARTED
                    sm_text_steps_properties.setAutoDraw(True)
                
                # if sm_text_steps_properties is active this frame...
                if sm_text_steps_properties.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_message* updates
                
                # if sm_text_message is starting this frame...
                if sm_text_message.status == NOT_STARTED and message_show:
                    # keep track of start time/frame for later
                    sm_text_message.frameNStart = frameN  # exact frame index
                    sm_text_message.tStart = t  # local t and not account for scr refresh
                    sm_text_message.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_message, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_message.status = STARTED
                    sm_text_message.setAutoDraw(True)
                
                # if sm_text_message is active this frame...
                if sm_text_message.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_instruction* updates
                
                # if sm_text_instruction is starting this frame...
                if sm_text_instruction.status == NOT_STARTED and instruction_show:
                    # keep track of start time/frame for later
                    sm_text_instruction.frameNStart = frameN  # exact frame index
                    sm_text_instruction.tStart = t  # local t and not account for scr refresh
                    sm_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_instruction, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_instruction.status = STARTED
                    sm_text_instruction.setAutoDraw(True)
                
                # if sm_text_instruction is active this frame...
                if sm_text_instruction.status == STARTED:
                    # update params
                    pass
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer], 
                        playbackComponents=[]
                    )
                    # skip the frame we paused on
                    continue
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    show_message.forceEnded = routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in show_message.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "show_message" ---
            for thisComponent in show_message.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for show_message
            show_message.tStop = globalClock.getTime(format='float')
            show_message.tStopRefresh = tThisFlipGlobal
            thisExp.addData('show_message.stopped', show_message.tStop)
            # Run 'End Routine' code from sm_code_update
            # Ick, but whatever. Makes sure that we don't ever
            # drop unexpectedly out of a move loop without 
            # resetting
            trial_start = True
            # the Routine "show_message" was not non-slip safe, so reset the non-slip timer
            routineTimer.reset()
        # completed 1.0 repeats of 'instructions_loop'
        
        
        # set up handler to look after randomisation of conditions etc
        online_trial_loop = data.TrialHandler2(
            name='online_trial_loop',
            nReps=online_trials_count, 
            method='fullRandom', 
            extraInfo=expInfo, 
            originPath=-1, 
            trialList=data.importConditions(snr_conditions_file), 
            seed=None, 
        )
        thisExp.addLoop(online_trial_loop)  # add the loop to the experiment
        thisOnline_trial_loop = online_trial_loop.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb = thisOnline_trial_loop.rgb)
        if thisOnline_trial_loop != None:
            for paramName in thisOnline_trial_loop:
                globals()[paramName] = thisOnline_trial_loop[paramName]
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        
        for thisOnline_trial_loop in online_trial_loop:
            currentLoop = online_trial_loop
            thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
            if thisSession is not None:
                # if running in a Session with a Liaison client, send data up to now
                thisSession.sendExperimentData()
            # abbreviate parameter names if possible (e.g. rgb = thisOnline_trial_loop.rgb)
            if thisOnline_trial_loop != None:
                for paramName in thisOnline_trial_loop:
                    globals()[paramName] = thisOnline_trial_loop[paramName]
            
            # set up handler to look after randomisation of conditions etc
            online_move_loop = data.TrialHandler2(
                name='online_move_loop',
                nReps=STEPS_MAX_PER_TRIAL, 
                method='sequential', 
                extraInfo=expInfo, 
                originPath=-1, 
                trialList=[None], 
                seed=None, 
            )
            thisExp.addLoop(online_move_loop)  # add the loop to the experiment
            thisOnline_move_loop = online_move_loop.trialList[0]  # so we can initialise stimuli with some values
            # abbreviate parameter names if possible (e.g. rgb = thisOnline_move_loop.rgb)
            if thisOnline_move_loop != None:
                for paramName in thisOnline_move_loop:
                    globals()[paramName] = thisOnline_move_loop[paramName]
            
            for thisOnline_move_loop in online_move_loop:
                currentLoop = online_move_loop
                thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
                # abbreviate parameter names if possible (e.g. rgb = thisOnline_move_loop.rgb)
                if thisOnline_move_loop != None:
                    for paramName in thisOnline_move_loop:
                        globals()[paramName] = thisOnline_move_loop[paramName]
                
                # --- Prepare to start Routine "online_task_run" ---
                # create an object to store info about Routine online_task_run
                online_task_run = data.Routine(
                    name='online_task_run',
                    components=[otr_kbd_response, otr_polygon_abscissa, otr_polygon_ordinate, otr_polygon_petey_midpoint, otr_image_pigeon, otr_image_right_seeds, otr_image_left_seeds, otr_polygon_coins_bar, otr_polygon_steps_bar, otr_polygon_status_bar_apex, otr_polygon_status_bar_base, otr_text_coins, otr_text_coins_label, otr_text_coins_properties, otr_text_steps, otr_text_steps_label, otr_text_steps_properties, otr_text_instruction],
                )
                online_task_run.status = NOT_STARTED
                continueRoutine = True
                # update component parameters for each repeat
                # create starting attributes for otr_kbd_response
                otr_kbd_response.keys = []
                otr_kbd_response.rt = []
                _otr_kbd_response_allKeys = []
                # Run 'Begin Routine' code from otr_code_update
                if trial_start is True:    
                    
                    # Clear keyboard events
                    otr_kbd_response.clearEvents()    
                    
                    # Default feedback
                    message_string = 'Failed to respond\nPlease respond more quickly'
                     
                    # Set instruction string
                    instruction_string = CHOICE_STRING
                    
                    # Initialize pigeon variables
                    pigeon_flip = 0.1
                    pigeon_true_x = ORIGIN_X
                    pigeon_shown_x = pigeon_true_x
                    pigeon_steps = [pigeon_shown_x]
                    
                    step_mean = ssc_step_mean
                    step_std = ssc_step_std
                
                    # Pay to play
                    coins_count += coins_paid_to_start_trial
                    steps_count += steps_taken_to_start_trial
                    
                    # Randomize direction, set std
                    if randint(0,2)==1:
                        correct_ans = 'right'
                        step_mean_val = abs(step_mean)
                    else:
                        correct_ans = 'left'
                        step_mean_val = -1*abs(step_mean)
                    step_std_val = step_std
                    
                    # Save trial parameters
                    thisExp.addData('step_mean_val', step_mean_val)
                    thisExp.addData('step_std_val', step_std_val)
                    thisExp.addData('block_number', block_number)
                    thisExp.addData('trial_number', trial_number)
                    thisExp.addData('correct', correct_ans)
                
                    # Set flag
                    trial_start = False
                    
                else:
                    
                    # Increment the counter and move the pigeon
                    steps_count += 1
                    move_pigeon(step_mean_val, step_std_val)
                    pigeon_steps.append(pigeon_shown_x)
                
                    # Update steps status
                    update_status(False, False, True)
                    
                # Check for end of task (no more steps)
                if steps_count >= steps_max: 
                        
                    # Drop out
                    endroutine = True;
                    online_move_loop.finished = True
                    online_trial_loop.finished = True
                otr_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
                otr_polygon_abscissa.setSize((EDGE_DISTANCE*2,LINE_HEIGHT))
                otr_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
                otr_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
                otr_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
                otr_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
                otr_image_pigeon.setPos((pigeon_shown_x,ORIGIN_Y))
                otr_image_pigeon.setSize((pigeon_flip, 0.1))
                otr_polygon_coins_bar.setPos((COINS_BAR_X,coins_bar_center))
                otr_polygon_coins_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
                otr_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
                otr_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
                otr_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
                otr_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
                otr_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
                otr_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
                otr_text_coins.setText(coins_count)
                otr_text_coins_properties.setText(coins_string)
                otr_text_steps.setText(steps_max - steps_count)
                otr_text_steps_properties.setText(steps_string)
                otr_text_instruction.setText(instruction_string)
                # store start times for online_task_run
                online_task_run.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
                online_task_run.tStart = globalClock.getTime(format='float')
                online_task_run.status = STARTED
                thisExp.addData('online_task_run.started', online_task_run.tStart)
                online_task_run.maxDuration = None
                # keep track of which components have finished
                online_task_runComponents = online_task_run.components
                for thisComponent in online_task_run.components:
                    thisComponent.tStart = None
                    thisComponent.tStop = None
                    thisComponent.tStartRefresh = None
                    thisComponent.tStopRefresh = None
                    if hasattr(thisComponent, 'status'):
                        thisComponent.status = NOT_STARTED
                # reset timers
                t = 0
                _timeToFirstFrame = win.getFutureFlipTime(clock="now")
                frameN = -1
                
                # --- Run Routine "online_task_run" ---
                # if trial has changed, end Routine now
                if isinstance(online_move_loop, data.TrialHandler2) and thisOnline_move_loop.thisN != online_move_loop.thisTrial.thisN:
                    continueRoutine = False
                online_task_run.forceEnded = routineForceEnded = not continueRoutine
                while continueRoutine:
                    # get current time
                    t = routineTimer.getTime()
                    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                    # update/draw components on each frame
                    
                    # *otr_kbd_response* updates
                    waitOnFlip = False
                    
                    # if otr_kbd_response is starting this frame...
                    if otr_kbd_response.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_kbd_response.frameNStart = frameN  # exact frame index
                        otr_kbd_response.tStart = t  # local t and not account for scr refresh
                        otr_kbd_response.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_kbd_response, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_kbd_response.status = STARTED
                        # keyboard checking is just starting
                        waitOnFlip = True
                        win.callOnFlip(otr_kbd_response.clock.reset)  # t=0 on next screen flip
                    
                    # if otr_kbd_response is stopping this frame...
                    if otr_kbd_response.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_kbd_response.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_kbd_response.tStop = t  # not accounting for scr refresh
                            otr_kbd_response.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_kbd_response.frameNStop = frameN  # exact frame index
                            # update status
                            otr_kbd_response.status = FINISHED
                            otr_kbd_response.status = FINISHED
                    if otr_kbd_response.status == STARTED and not waitOnFlip:
                        theseKeys = otr_kbd_response.getKeys(keyList=['left','right'], ignoreKeys=["escape"], waitRelease=False)
                        _otr_kbd_response_allKeys.extend(theseKeys)
                        if len(_otr_kbd_response_allKeys):
                            otr_kbd_response.keys = _otr_kbd_response_allKeys[-1].name  # just the last key pressed
                            otr_kbd_response.rt = _otr_kbd_response_allKeys[-1].rt
                            otr_kbd_response.duration = _otr_kbd_response_allKeys[-1].duration
                            # a response ends the routine
                            continueRoutine = False
                    
                    # *otr_polygon_abscissa* updates
                    
                    # if otr_polygon_abscissa is starting this frame...
                    if otr_polygon_abscissa.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_polygon_abscissa.frameNStart = frameN  # exact frame index
                        otr_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                        otr_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_polygon_abscissa.status = STARTED
                        otr_polygon_abscissa.setAutoDraw(True)
                    
                    # if otr_polygon_abscissa is active this frame...
                    if otr_polygon_abscissa.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_polygon_abscissa is stopping this frame...
                    if otr_polygon_abscissa.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_polygon_abscissa.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_polygon_abscissa.tStop = t  # not accounting for scr refresh
                            otr_polygon_abscissa.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_polygon_abscissa.frameNStop = frameN  # exact frame index
                            # update status
                            otr_polygon_abscissa.status = FINISHED
                            otr_polygon_abscissa.setAutoDraw(False)
                    
                    # *otr_polygon_ordinate* updates
                    
                    # if otr_polygon_ordinate is starting this frame...
                    if otr_polygon_ordinate.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_polygon_ordinate.frameNStart = frameN  # exact frame index
                        otr_polygon_ordinate.tStart = t  # local t and not account for scr refresh
                        otr_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_polygon_ordinate.status = STARTED
                        otr_polygon_ordinate.setAutoDraw(True)
                    
                    # if otr_polygon_ordinate is active this frame...
                    if otr_polygon_ordinate.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_polygon_ordinate is stopping this frame...
                    if otr_polygon_ordinate.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_polygon_ordinate.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_polygon_ordinate.tStop = t  # not accounting for scr refresh
                            otr_polygon_ordinate.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_polygon_ordinate.frameNStop = frameN  # exact frame index
                            # update status
                            otr_polygon_ordinate.status = FINISHED
                            otr_polygon_ordinate.setAutoDraw(False)
                    
                    # *otr_polygon_petey_midpoint* updates
                    
                    # if otr_polygon_petey_midpoint is starting this frame...
                    if otr_polygon_petey_midpoint.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
                        otr_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
                        otr_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_polygon_petey_midpoint.status = STARTED
                        otr_polygon_petey_midpoint.setAutoDraw(True)
                    
                    # if otr_polygon_petey_midpoint is active this frame...
                    if otr_polygon_petey_midpoint.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_polygon_petey_midpoint is stopping this frame...
                    if otr_polygon_petey_midpoint.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_polygon_petey_midpoint.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_polygon_petey_midpoint.tStop = t  # not accounting for scr refresh
                            otr_polygon_petey_midpoint.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_polygon_petey_midpoint.frameNStop = frameN  # exact frame index
                            # update status
                            otr_polygon_petey_midpoint.status = FINISHED
                            otr_polygon_petey_midpoint.setAutoDraw(False)
                    
                    # *otr_image_pigeon* updates
                    
                    # if otr_image_pigeon is starting this frame...
                    if otr_image_pigeon.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_image_pigeon.frameNStart = frameN  # exact frame index
                        otr_image_pigeon.tStart = t  # local t and not account for scr refresh
                        otr_image_pigeon.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_image_pigeon, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_image_pigeon.status = STARTED
                        otr_image_pigeon.setAutoDraw(True)
                    
                    # if otr_image_pigeon is active this frame...
                    if otr_image_pigeon.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_image_pigeon is stopping this frame...
                    if otr_image_pigeon.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_image_pigeon.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_image_pigeon.tStop = t  # not accounting for scr refresh
                            otr_image_pigeon.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_image_pigeon.frameNStop = frameN  # exact frame index
                            # update status
                            otr_image_pigeon.status = FINISHED
                            otr_image_pigeon.setAutoDraw(False)
                    
                    # *otr_image_right_seeds* updates
                    
                    # if otr_image_right_seeds is starting this frame...
                    if otr_image_right_seeds.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_image_right_seeds.frameNStart = frameN  # exact frame index
                        otr_image_right_seeds.tStart = t  # local t and not account for scr refresh
                        otr_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_image_right_seeds.status = STARTED
                        otr_image_right_seeds.setAutoDraw(True)
                    
                    # if otr_image_right_seeds is active this frame...
                    if otr_image_right_seeds.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_image_right_seeds is stopping this frame...
                    if otr_image_right_seeds.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_image_right_seeds.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_image_right_seeds.tStop = t  # not accounting for scr refresh
                            otr_image_right_seeds.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_image_right_seeds.frameNStop = frameN  # exact frame index
                            # update status
                            otr_image_right_seeds.status = FINISHED
                            otr_image_right_seeds.setAutoDraw(False)
                    
                    # *otr_image_left_seeds* updates
                    
                    # if otr_image_left_seeds is starting this frame...
                    if otr_image_left_seeds.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_image_left_seeds.frameNStart = frameN  # exact frame index
                        otr_image_left_seeds.tStart = t  # local t and not account for scr refresh
                        otr_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_image_left_seeds.status = STARTED
                        otr_image_left_seeds.setAutoDraw(True)
                    
                    # if otr_image_left_seeds is active this frame...
                    if otr_image_left_seeds.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_image_left_seeds is stopping this frame...
                    if otr_image_left_seeds.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_image_left_seeds.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_image_left_seeds.tStop = t  # not accounting for scr refresh
                            otr_image_left_seeds.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_image_left_seeds.frameNStop = frameN  # exact frame index
                            # update status
                            otr_image_left_seeds.status = FINISHED
                            otr_image_left_seeds.setAutoDraw(False)
                    
                    # *otr_polygon_coins_bar* updates
                    
                    # if otr_polygon_coins_bar is starting this frame...
                    if otr_polygon_coins_bar.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_polygon_coins_bar.frameNStart = frameN  # exact frame index
                        otr_polygon_coins_bar.tStart = t  # local t and not account for scr refresh
                        otr_polygon_coins_bar.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_polygon_coins_bar, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_polygon_coins_bar.status = STARTED
                        otr_polygon_coins_bar.setAutoDraw(True)
                    
                    # if otr_polygon_coins_bar is active this frame...
                    if otr_polygon_coins_bar.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_polygon_coins_bar is stopping this frame...
                    if otr_polygon_coins_bar.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_polygon_coins_bar.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_polygon_coins_bar.tStop = t  # not accounting for scr refresh
                            otr_polygon_coins_bar.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_polygon_coins_bar.frameNStop = frameN  # exact frame index
                            # update status
                            otr_polygon_coins_bar.status = FINISHED
                            otr_polygon_coins_bar.setAutoDraw(False)
                    
                    # *otr_polygon_steps_bar* updates
                    
                    # if otr_polygon_steps_bar is starting this frame...
                    if otr_polygon_steps_bar.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_polygon_steps_bar.frameNStart = frameN  # exact frame index
                        otr_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
                        otr_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_polygon_steps_bar.status = STARTED
                        otr_polygon_steps_bar.setAutoDraw(True)
                    
                    # if otr_polygon_steps_bar is active this frame...
                    if otr_polygon_steps_bar.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_polygon_steps_bar is stopping this frame...
                    if otr_polygon_steps_bar.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_polygon_steps_bar.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_polygon_steps_bar.tStop = t  # not accounting for scr refresh
                            otr_polygon_steps_bar.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_polygon_steps_bar.frameNStop = frameN  # exact frame index
                            # update status
                            otr_polygon_steps_bar.status = FINISHED
                            otr_polygon_steps_bar.setAutoDraw(False)
                    
                    # *otr_polygon_status_bar_apex* updates
                    
                    # if otr_polygon_status_bar_apex is starting this frame...
                    if otr_polygon_status_bar_apex.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
                        otr_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
                        otr_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_polygon_status_bar_apex.status = STARTED
                        otr_polygon_status_bar_apex.setAutoDraw(True)
                    
                    # if otr_polygon_status_bar_apex is active this frame...
                    if otr_polygon_status_bar_apex.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_polygon_status_bar_apex is stopping this frame...
                    if otr_polygon_status_bar_apex.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_polygon_status_bar_apex.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_polygon_status_bar_apex.tStop = t  # not accounting for scr refresh
                            otr_polygon_status_bar_apex.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_polygon_status_bar_apex.frameNStop = frameN  # exact frame index
                            # update status
                            otr_polygon_status_bar_apex.status = FINISHED
                            otr_polygon_status_bar_apex.setAutoDraw(False)
                    
                    # *otr_polygon_status_bar_base* updates
                    
                    # if otr_polygon_status_bar_base is starting this frame...
                    if otr_polygon_status_bar_base.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_polygon_status_bar_base.frameNStart = frameN  # exact frame index
                        otr_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
                        otr_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_polygon_status_bar_base.status = STARTED
                        otr_polygon_status_bar_base.setAutoDraw(True)
                    
                    # if otr_polygon_status_bar_base is active this frame...
                    if otr_polygon_status_bar_base.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_polygon_status_bar_base is stopping this frame...
                    if otr_polygon_status_bar_base.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_polygon_status_bar_base.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_polygon_status_bar_base.tStop = t  # not accounting for scr refresh
                            otr_polygon_status_bar_base.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_polygon_status_bar_base.frameNStop = frameN  # exact frame index
                            # update status
                            otr_polygon_status_bar_base.status = FINISHED
                            otr_polygon_status_bar_base.setAutoDraw(False)
                    
                    # *otr_text_coins* updates
                    
                    # if otr_text_coins is starting this frame...
                    if otr_text_coins.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_text_coins.frameNStart = frameN  # exact frame index
                        otr_text_coins.tStart = t  # local t and not account for scr refresh
                        otr_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_text_coins, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_text_coins.status = STARTED
                        otr_text_coins.setAutoDraw(True)
                    
                    # if otr_text_coins is active this frame...
                    if otr_text_coins.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_text_coins is stopping this frame...
                    if otr_text_coins.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_text_coins.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_text_coins.tStop = t  # not accounting for scr refresh
                            otr_text_coins.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_text_coins.frameNStop = frameN  # exact frame index
                            # update status
                            otr_text_coins.status = FINISHED
                            otr_text_coins.setAutoDraw(False)
                    
                    # *otr_text_coins_label* updates
                    
                    # if otr_text_coins_label is starting this frame...
                    if otr_text_coins_label.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_text_coins_label.frameNStart = frameN  # exact frame index
                        otr_text_coins_label.tStart = t  # local t and not account for scr refresh
                        otr_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_text_coins_label, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_text_coins_label.status = STARTED
                        otr_text_coins_label.setAutoDraw(True)
                    
                    # if otr_text_coins_label is active this frame...
                    if otr_text_coins_label.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_text_coins_label is stopping this frame...
                    if otr_text_coins_label.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_text_coins_label.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_text_coins_label.tStop = t  # not accounting for scr refresh
                            otr_text_coins_label.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_text_coins_label.frameNStop = frameN  # exact frame index
                            # update status
                            otr_text_coins_label.status = FINISHED
                            otr_text_coins_label.setAutoDraw(False)
                    
                    # *otr_text_coins_properties* updates
                    
                    # if otr_text_coins_properties is starting this frame...
                    if otr_text_coins_properties.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_text_coins_properties.frameNStart = frameN  # exact frame index
                        otr_text_coins_properties.tStart = t  # local t and not account for scr refresh
                        otr_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_text_coins_properties.status = STARTED
                        otr_text_coins_properties.setAutoDraw(True)
                    
                    # if otr_text_coins_properties is active this frame...
                    if otr_text_coins_properties.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_text_coins_properties is stopping this frame...
                    if otr_text_coins_properties.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_text_coins_properties.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_text_coins_properties.tStop = t  # not accounting for scr refresh
                            otr_text_coins_properties.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_text_coins_properties.frameNStop = frameN  # exact frame index
                            # update status
                            otr_text_coins_properties.status = FINISHED
                            otr_text_coins_properties.setAutoDraw(False)
                    
                    # *otr_text_steps* updates
                    
                    # if otr_text_steps is starting this frame...
                    if otr_text_steps.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_text_steps.frameNStart = frameN  # exact frame index
                        otr_text_steps.tStart = t  # local t and not account for scr refresh
                        otr_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_text_steps, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_text_steps.status = STARTED
                        otr_text_steps.setAutoDraw(True)
                    
                    # if otr_text_steps is active this frame...
                    if otr_text_steps.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_text_steps is stopping this frame...
                    if otr_text_steps.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_text_steps.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_text_steps.tStop = t  # not accounting for scr refresh
                            otr_text_steps.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_text_steps.frameNStop = frameN  # exact frame index
                            # update status
                            otr_text_steps.status = FINISHED
                            otr_text_steps.setAutoDraw(False)
                    
                    # *otr_text_steps_label* updates
                    
                    # if otr_text_steps_label is starting this frame...
                    if otr_text_steps_label.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_text_steps_label.frameNStart = frameN  # exact frame index
                        otr_text_steps_label.tStart = t  # local t and not account for scr refresh
                        otr_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_text_steps_label, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_text_steps_label.status = STARTED
                        otr_text_steps_label.setAutoDraw(True)
                    
                    # if otr_text_steps_label is active this frame...
                    if otr_text_steps_label.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_text_steps_label is stopping this frame...
                    if otr_text_steps_label.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_text_steps_label.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_text_steps_label.tStop = t  # not accounting for scr refresh
                            otr_text_steps_label.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_text_steps_label.frameNStop = frameN  # exact frame index
                            # update status
                            otr_text_steps_label.status = FINISHED
                            otr_text_steps_label.setAutoDraw(False)
                    
                    # *otr_text_steps_properties* updates
                    
                    # if otr_text_steps_properties is starting this frame...
                    if otr_text_steps_properties.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_text_steps_properties.frameNStart = frameN  # exact frame index
                        otr_text_steps_properties.tStart = t  # local t and not account for scr refresh
                        otr_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_text_steps_properties.status = STARTED
                        otr_text_steps_properties.setAutoDraw(True)
                    
                    # if otr_text_steps_properties is active this frame...
                    if otr_text_steps_properties.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_text_steps_properties is stopping this frame...
                    if otr_text_steps_properties.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_text_steps_properties.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_text_steps_properties.tStop = t  # not accounting for scr refresh
                            otr_text_steps_properties.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_text_steps_properties.frameNStop = frameN  # exact frame index
                            # update status
                            otr_text_steps_properties.status = FINISHED
                            otr_text_steps_properties.setAutoDraw(False)
                    
                    # *otr_text_instruction* updates
                    
                    # if otr_text_instruction is starting this frame...
                    if otr_text_instruction.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        otr_text_instruction.frameNStart = frameN  # exact frame index
                        otr_text_instruction.tStart = t  # local t and not account for scr refresh
                        otr_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(otr_text_instruction, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        otr_text_instruction.status = STARTED
                        otr_text_instruction.setAutoDraw(True)
                    
                    # if otr_text_instruction is active this frame...
                    if otr_text_instruction.status == STARTED:
                        # update params
                        pass
                    
                    # if otr_text_instruction is stopping this frame...
                    if otr_text_instruction.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > otr_text_instruction.tStartRefresh + UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            otr_text_instruction.tStop = t  # not accounting for scr refresh
                            otr_text_instruction.tStopRefresh = tThisFlipGlobal  # on global time
                            otr_text_instruction.frameNStop = frameN  # exact frame index
                            # update status
                            otr_text_instruction.status = FINISHED
                            otr_text_instruction.setAutoDraw(False)
                    
                    # check for quit (typically the Esc key)
                    if defaultKeyboard.getKeys(keyList=["escape"]):
                        thisExp.status = FINISHED
                    if thisExp.status == FINISHED or endExpNow:
                        endExperiment(thisExp, win=win)
                        return
                    # pause experiment here if requested
                    if thisExp.status == PAUSED:
                        pauseExperiment(
                            thisExp=thisExp, 
                            win=win, 
                            timers=[routineTimer], 
                            playbackComponents=[]
                        )
                        # skip the frame we paused on
                        continue
                    
                    # check if all components have finished
                    if not continueRoutine:  # a component has requested a forced-end of Routine
                        online_task_run.forceEnded = routineForceEnded = True
                        break
                    continueRoutine = False  # will revert to True if at least one component still running
                    for thisComponent in online_task_run.components:
                        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                            continueRoutine = True
                            break  # at least one component has not yet finished
                    
                    # refresh the screen
                    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                        win.flip()
                
                # --- Ending Routine "online_task_run" ---
                for thisComponent in online_task_run.components:
                    if hasattr(thisComponent, "setAutoDraw"):
                        thisComponent.setAutoDraw(False)
                # store stop times for online_task_run
                online_task_run.tStop = globalClock.getTime(format='float')
                online_task_run.tStopRefresh = tThisFlipGlobal
                thisExp.addData('online_task_run.stopped', online_task_run.tStop)
                # Run 'End Routine' code from otr_code_update
                if online_trial_loop.finished is True:
                    
                     # Cap # of steps
                    steps_count = steps_max;
                    
                    # Save results
                    thisExp.addData('pigeon_steps', pigeon_steps)
                    thisExp.addData('coins_count', coins_count)
                    thisExp.addData('steps_count', steps_count)
                    thisExp.addData('choice', 'none')
                
                    # Drop out of this block loop
                    set_block_end_string()
                    thisExp.addData('bonus', bonus_count)
                
                elif otr_kbd_response.keys:
                        
                    # Parse answer and update status
                    if correct_ans in otr_kbd_response.keys:
                        coins_count += coins_gained_per_correct
                        set_trial_end_string(CORRECT_STRING)
                        update_status(False, True, False)
                    else:
                        coins_count -= coins_lost_per_error
                        steps_count += steps_lost_per_error
                        if steps_count >= steps_max:
                            steps_count = steps_max
                        set_trial_end_string(ERROR_STRING)
                    
                    # Save results
                    thisExp.addData('pigeon_steps', pigeon_steps)
                    thisExp.addData('coins_count', coins_count)
                    thisExp.addData('steps_count', steps_count)    
                    thisExp.addData('choice', otr_kbd_response.keys)
                    
                    # update status
                    update_status(False, True, True)
                            
                    # Clear keyboard events
                    otr_kbd_response.clearEvents()    
                
                    # Drop out of pigeon movement loop
                    trial_number += 1
                    online_move_loop.finished = True
                    
                    
                
                # the Routine "online_task_run" was not non-slip safe, so reset the non-slip timer
                routineTimer.reset()
            # completed STEPS_MAX_PER_TRIAL repeats of 'online_move_loop'
            
            
            # --- Prepare to start Routine "show_message" ---
            # create an object to store info about Routine show_message
            show_message = data.Routine(
                name='show_message',
                components=[sm_kbd_continue, sm_polygon_box, sm_polygon_abscissa, sm_polygon_ordinate, sm_polygon_petey_midpoint, sm_image_right_seeds, sm_image_left_seeds, sm_image_petey, sm_polygon_coin_bar, sm_polygon_steps_bar, sm_polygon_status_bar_apex, sm_polygon_status_bar_base, sm_text_coins, sm_text_coins_label, sm_text_coins_properties, sm_text_steps, sm_text_steps_label, sm_text_steps_properties, sm_text_message, sm_text_instruction],
            )
            show_message.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            # create starting attributes for sm_kbd_continue
            sm_kbd_continue.keys = []
            sm_kbd_continue.rt = []
            _sm_kbd_continue_allKeys = []
            # Run 'Begin Routine' code from sm_code_update
            # Clear all keyboard events
            sm_kbd_continue.clearEvents()
            sm_polygon_box.setPos((box_x,box_y))
            sm_polygon_box.setSize((box_width,box_height))
            sm_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
            sm_polygon_abscissa.setSize((ABSCISSA_WIDTH,LINE_HEIGHT))
            sm_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
            sm_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
            sm_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
            sm_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
            sm_image_petey.setPos((pigeon_shown_x,ORIGIN_Y))
            sm_image_petey.setSize((pigeon_flip, 0.1))
            sm_polygon_coin_bar.setPos((COINS_BAR_X,coins_bar_center))
            sm_polygon_coin_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
            sm_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
            sm_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
            sm_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
            sm_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
            sm_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
            sm_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
            sm_text_coins.setText(coins_count)
            sm_text_coins_properties.setText(coins_string)
            sm_text_steps.setText(steps_max - steps_count)
            sm_text_steps_properties.setText(steps_string)
            sm_text_message.setText(message_string)
            sm_text_instruction.setText(CONTINUE_STRING)
            # store start times for show_message
            show_message.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            show_message.tStart = globalClock.getTime(format='float')
            show_message.status = STARTED
            thisExp.addData('show_message.started', show_message.tStart)
            show_message.maxDuration = None
            # keep track of which components have finished
            show_messageComponents = show_message.components
            for thisComponent in show_message.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "show_message" ---
            # if trial has changed, end Routine now
            if isinstance(online_trial_loop, data.TrialHandler2) and thisOnline_trial_loop.thisN != online_trial_loop.thisTrial.thisN:
                continueRoutine = False
            show_message.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                # *sm_kbd_continue* updates
                waitOnFlip = False
                
                # if sm_kbd_continue is starting this frame...
                if sm_kbd_continue.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    sm_kbd_continue.frameNStart = frameN  # exact frame index
                    sm_kbd_continue.tStart = t  # local t and not account for scr refresh
                    sm_kbd_continue.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_kbd_continue, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_kbd_continue.status = STARTED
                    # keyboard checking is just starting
                    waitOnFlip = True
                    win.callOnFlip(sm_kbd_continue.clock.reset)  # t=0 on next screen flip
                    win.callOnFlip(sm_kbd_continue.clearEvents, eventType='keyboard')  # clear events on next screen flip
                if sm_kbd_continue.status == STARTED and not waitOnFlip:
                    theseKeys = sm_kbd_continue.getKeys(keyList=['space', 'up', 'down', 'left', 'right'], ignoreKeys=["escape"], waitRelease=False)
                    _sm_kbd_continue_allKeys.extend(theseKeys)
                    if len(_sm_kbd_continue_allKeys):
                        sm_kbd_continue.keys = _sm_kbd_continue_allKeys[-1].name  # just the last key pressed
                        sm_kbd_continue.rt = _sm_kbd_continue_allKeys[-1].rt
                        sm_kbd_continue.duration = _sm_kbd_continue_allKeys[-1].duration
                        # a response ends the routine
                        continueRoutine = False
                
                # *sm_polygon_box* updates
                
                # if sm_polygon_box is starting this frame...
                if sm_polygon_box.status == NOT_STARTED and box_show:
                    # keep track of start time/frame for later
                    sm_polygon_box.frameNStart = frameN  # exact frame index
                    sm_polygon_box.tStart = t  # local t and not account for scr refresh
                    sm_polygon_box.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_box, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_box.status = STARTED
                    sm_polygon_box.setAutoDraw(True)
                
                # if sm_polygon_box is active this frame...
                if sm_polygon_box.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_abscissa* updates
                
                # if sm_polygon_abscissa is starting this frame...
                if sm_polygon_abscissa.status == NOT_STARTED and pigeon_show:
                    # keep track of start time/frame for later
                    sm_polygon_abscissa.frameNStart = frameN  # exact frame index
                    sm_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                    sm_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_abscissa.status = STARTED
                    sm_polygon_abscissa.setAutoDraw(True)
                
                # if sm_polygon_abscissa is active this frame...
                if sm_polygon_abscissa.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_ordinate* updates
                
                # if sm_polygon_ordinate is starting this frame...
                if sm_polygon_ordinate.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_ordinate.frameNStart = frameN  # exact frame index
                    sm_polygon_ordinate.tStart = t  # local t and not account for scr refresh
                    sm_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_ordinate.status = STARTED
                    sm_polygon_ordinate.setAutoDraw(True)
                
                # if sm_polygon_ordinate is active this frame...
                if sm_polygon_ordinate.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_petey_midpoint* updates
                
                # if sm_polygon_petey_midpoint is starting this frame...
                if sm_polygon_petey_midpoint.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
                    sm_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
                    sm_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_petey_midpoint.status = STARTED
                    sm_polygon_petey_midpoint.setAutoDraw(True)
                
                # if sm_polygon_petey_midpoint is active this frame...
                if sm_polygon_petey_midpoint.status == STARTED:
                    # update params
                    pass
                
                # *sm_image_right_seeds* updates
                
                # if sm_image_right_seeds is starting this frame...
                if sm_image_right_seeds.status == NOT_STARTED and seeds_show:
                    # keep track of start time/frame for later
                    sm_image_right_seeds.frameNStart = frameN  # exact frame index
                    sm_image_right_seeds.tStart = t  # local t and not account for scr refresh
                    sm_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_image_right_seeds.status = STARTED
                    sm_image_right_seeds.setAutoDraw(True)
                
                # if sm_image_right_seeds is active this frame...
                if sm_image_right_seeds.status == STARTED:
                    # update params
                    pass
                
                # *sm_image_left_seeds* updates
                
                # if sm_image_left_seeds is starting this frame...
                if sm_image_left_seeds.status == NOT_STARTED and seeds_show:
                    # keep track of start time/frame for later
                    sm_image_left_seeds.frameNStart = frameN  # exact frame index
                    sm_image_left_seeds.tStart = t  # local t and not account for scr refresh
                    sm_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_image_left_seeds.status = STARTED
                    sm_image_left_seeds.setAutoDraw(True)
                
                # if sm_image_left_seeds is active this frame...
                if sm_image_left_seeds.status == STARTED:
                    # update params
                    pass
                
                # *sm_image_petey* updates
                
                # if sm_image_petey is starting this frame...
                if sm_image_petey.status == NOT_STARTED and pigeon_show:
                    # keep track of start time/frame for later
                    sm_image_petey.frameNStart = frameN  # exact frame index
                    sm_image_petey.tStart = t  # local t and not account for scr refresh
                    sm_image_petey.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_image_petey, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_image_petey.status = STARTED
                    sm_image_petey.setAutoDraw(True)
                
                # if sm_image_petey is active this frame...
                if sm_image_petey.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_coin_bar* updates
                
                # if sm_polygon_coin_bar is starting this frame...
                if sm_polygon_coin_bar.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_coin_bar.frameNStart = frameN  # exact frame index
                    sm_polygon_coin_bar.tStart = t  # local t and not account for scr refresh
                    sm_polygon_coin_bar.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_coin_bar, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_coin_bar.status = STARTED
                    sm_polygon_coin_bar.setAutoDraw(True)
                
                # if sm_polygon_coin_bar is active this frame...
                if sm_polygon_coin_bar.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_steps_bar* updates
                
                # if sm_polygon_steps_bar is starting this frame...
                if sm_polygon_steps_bar.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_steps_bar.frameNStart = frameN  # exact frame index
                    sm_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
                    sm_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_steps_bar.status = STARTED
                    sm_polygon_steps_bar.setAutoDraw(True)
                
                # if sm_polygon_steps_bar is active this frame...
                if sm_polygon_steps_bar.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_status_bar_apex* updates
                
                # if sm_polygon_status_bar_apex is starting this frame...
                if sm_polygon_status_bar_apex.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
                    sm_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
                    sm_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_status_bar_apex.status = STARTED
                    sm_polygon_status_bar_apex.setAutoDraw(True)
                
                # if sm_polygon_status_bar_apex is active this frame...
                if sm_polygon_status_bar_apex.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_status_bar_base* updates
                
                # if sm_polygon_status_bar_base is starting this frame...
                if sm_polygon_status_bar_base.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_status_bar_base.frameNStart = frameN  # exact frame index
                    sm_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
                    sm_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_status_bar_base.status = STARTED
                    sm_polygon_status_bar_base.setAutoDraw(True)
                
                # if sm_polygon_status_bar_base is active this frame...
                if sm_polygon_status_bar_base.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_coins* updates
                
                # if sm_text_coins is starting this frame...
                if sm_text_coins.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_coins.frameNStart = frameN  # exact frame index
                    sm_text_coins.tStart = t  # local t and not account for scr refresh
                    sm_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_coins, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_coins.status = STARTED
                    sm_text_coins.setAutoDraw(True)
                
                # if sm_text_coins is active this frame...
                if sm_text_coins.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_coins_label* updates
                
                # if sm_text_coins_label is starting this frame...
                if sm_text_coins_label.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_coins_label.frameNStart = frameN  # exact frame index
                    sm_text_coins_label.tStart = t  # local t and not account for scr refresh
                    sm_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_coins_label, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_coins_label.status = STARTED
                    sm_text_coins_label.setAutoDraw(True)
                
                # if sm_text_coins_label is active this frame...
                if sm_text_coins_label.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_coins_properties* updates
                
                # if sm_text_coins_properties is starting this frame...
                if sm_text_coins_properties.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_coins_properties.frameNStart = frameN  # exact frame index
                    sm_text_coins_properties.tStart = t  # local t and not account for scr refresh
                    sm_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_coins_properties.status = STARTED
                    sm_text_coins_properties.setAutoDraw(True)
                
                # if sm_text_coins_properties is active this frame...
                if sm_text_coins_properties.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_steps* updates
                
                # if sm_text_steps is starting this frame...
                if sm_text_steps.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_steps.frameNStart = frameN  # exact frame index
                    sm_text_steps.tStart = t  # local t and not account for scr refresh
                    sm_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_steps, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_steps.status = STARTED
                    sm_text_steps.setAutoDraw(True)
                
                # if sm_text_steps is active this frame...
                if sm_text_steps.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_steps_label* updates
                
                # if sm_text_steps_label is starting this frame...
                if sm_text_steps_label.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_steps_label.frameNStart = frameN  # exact frame index
                    sm_text_steps_label.tStart = t  # local t and not account for scr refresh
                    sm_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_steps_label, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_steps_label.status = STARTED
                    sm_text_steps_label.setAutoDraw(True)
                
                # if sm_text_steps_label is active this frame...
                if sm_text_steps_label.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_steps_properties* updates
                
                # if sm_text_steps_properties is starting this frame...
                if sm_text_steps_properties.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_steps_properties.frameNStart = frameN  # exact frame index
                    sm_text_steps_properties.tStart = t  # local t and not account for scr refresh
                    sm_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_steps_properties.status = STARTED
                    sm_text_steps_properties.setAutoDraw(True)
                
                # if sm_text_steps_properties is active this frame...
                if sm_text_steps_properties.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_message* updates
                
                # if sm_text_message is starting this frame...
                if sm_text_message.status == NOT_STARTED and message_show:
                    # keep track of start time/frame for later
                    sm_text_message.frameNStart = frameN  # exact frame index
                    sm_text_message.tStart = t  # local t and not account for scr refresh
                    sm_text_message.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_message, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_message.status = STARTED
                    sm_text_message.setAutoDraw(True)
                
                # if sm_text_message is active this frame...
                if sm_text_message.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_instruction* updates
                
                # if sm_text_instruction is starting this frame...
                if sm_text_instruction.status == NOT_STARTED and instruction_show:
                    # keep track of start time/frame for later
                    sm_text_instruction.frameNStart = frameN  # exact frame index
                    sm_text_instruction.tStart = t  # local t and not account for scr refresh
                    sm_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_instruction, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_instruction.status = STARTED
                    sm_text_instruction.setAutoDraw(True)
                
                # if sm_text_instruction is active this frame...
                if sm_text_instruction.status == STARTED:
                    # update params
                    pass
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer], 
                        playbackComponents=[]
                    )
                    # skip the frame we paused on
                    continue
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    show_message.forceEnded = routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in show_message.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "show_message" ---
            for thisComponent in show_message.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for show_message
            show_message.tStop = globalClock.getTime(format='float')
            show_message.tStopRefresh = tThisFlipGlobal
            thisExp.addData('show_message.stopped', show_message.tStop)
            # Run 'End Routine' code from sm_code_update
            # Ick, but whatever. Makes sure that we don't ever
            # drop unexpectedly out of a move loop without 
            # resetting
            trial_start = True
            # the Routine "show_message" was not non-slip safe, so reset the non-slip timer
            routineTimer.reset()
            thisExp.nextEntry()
            
        # completed online_trials_count repeats of 'online_trial_loop'
        
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        
        # set up handler to look after randomisation of conditions etc
        predefined_trial_loop = data.TrialHandler2(
            name='predefined_trial_loop',
            nReps=predefined_trials_count, 
            method='random', 
            extraInfo=expInfo, 
            originPath=-1, 
            trialList=data.importConditions(snr_conditions_file), 
            seed=None, 
        )
        thisExp.addLoop(predefined_trial_loop)  # add the loop to the experiment
        thisPredefined_trial_loop = predefined_trial_loop.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb = thisPredefined_trial_loop.rgb)
        if thisPredefined_trial_loop != None:
            for paramName in thisPredefined_trial_loop:
                globals()[paramName] = thisPredefined_trial_loop[paramName]
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        
        for thisPredefined_trial_loop in predefined_trial_loop:
            currentLoop = predefined_trial_loop
            thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
            if thisSession is not None:
                # if running in a Session with a Liaison client, send data up to now
                thisSession.sendExperimentData()
            # abbreviate parameter names if possible (e.g. rgb = thisPredefined_trial_loop.rgb)
            if thisPredefined_trial_loop != None:
                for paramName in thisPredefined_trial_loop:
                    globals()[paramName] = thisPredefined_trial_loop[paramName]
            
            # set up handler to look after randomisation of conditions etc
            predefined_setup_loop = data.TrialHandler2(
                name='predefined_setup_loop',
                nReps=1000.0, 
                method='sequential', 
                extraInfo=expInfo, 
                originPath=-1, 
                trialList=[None], 
                seed=None, 
            )
            thisExp.addLoop(predefined_setup_loop)  # add the loop to the experiment
            thisPredefined_setup_loop = predefined_setup_loop.trialList[0]  # so we can initialise stimuli with some values
            # abbreviate parameter names if possible (e.g. rgb = thisPredefined_setup_loop.rgb)
            if thisPredefined_setup_loop != None:
                for paramName in thisPredefined_setup_loop:
                    globals()[paramName] = thisPredefined_setup_loop[paramName]
            
            for thisPredefined_setup_loop in predefined_setup_loop:
                currentLoop = predefined_setup_loop
                thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
                # abbreviate parameter names if possible (e.g. rgb = thisPredefined_setup_loop.rgb)
                if thisPredefined_setup_loop != None:
                    for paramName in thisPredefined_setup_loop:
                        globals()[paramName] = thisPredefined_setup_loop[paramName]
                
                # --- Prepare to start Routine "predefined_task_setup" ---
                # create an object to store info about Routine predefined_task_setup
                predefined_task_setup = data.Routine(
                    name='predefined_task_setup',
                    components=[pts_kbd_update, pts_polygon_abscissa, pts_polygon_ordinate, pts_polygon_petey_midpoint, pts_polygon_left_bound, pts_polygon_right_bound, pts_polygon_stay_left, pts_polygon_stay_right, pts_image_pigeon, pts_image_right_seeds, pts_image_left_seeds, pts_polygon_coins_bar, pts_polygon_steps_bar, pts_polygon_status_bar_apex, pts_polygon_status_bar_base, pts_text_coins, pts_text_coins_label, pts_text_coins_properties, pts_text_steps, pts_text_steps_label, pts_text_steps_properties, pts_text_instruction, pts_text_hint, pts_text_hint_2],
                )
                predefined_task_setup.status = NOT_STARTED
                continueRoutine = True
                # update component parameters for each repeat
                # create starting attributes for pts_kbd_update
                pts_kbd_update.keys = []
                pts_kbd_update.rt = []
                _pts_kbd_update_allKeys = []
                # Run 'Begin Routine' code from pts_code_update
                if trial_start:
                
                    # Do this once
                    trial_start = False
                
                    # clear keyboard events
                    pts_kbd_update.clearEvents()    
                    
                    # Set instruction string
                    if task_type == 'predefined':
                        instruction_string = PREDEFINED_CHOICE_STRING
                    else:
                        if steps_count < min_steps_to_commit:
                            instruction_string = PREDEFINED_CHOICE_STRING
                        else:
                            instruction_string = PREDEFINED_CHOICE_COMMIT_STRING
                        
                    # Default feedback
                    message_string = PREDEFINED_NO_CHOICE_STRING
                    choice = 'max'
                
                    # initialize pigeon variables
                    pigeon_flip = 0.1
                    pigeon_true_x = 0
                    pigeon_shown_x = pigeon_true_x
                    pigeon_steps = [pigeon_shown_x]
                    predefined_bound_x_previous = predefined_bound_x
                
                    # Pay to play
                    coins_count += coins_paid_to_start_trial
                    steps_count += steps_taken_to_start_trial
                    
                    step_mean = ssc_step_mean
                    step_std = ssc_step_std
                    
                    # Randomize direction, save std
                    if randint(0,2)==1:
                        correct_ans = 'right'
                        step_mean_val = abs(step_mean)
                    else:
                        correct_ans = 'left'
                        step_mean_val = -1*abs(step_mean)
                    step_std_val = step_std
                
                    # Save trial parameters
                    thisExp.addData('step_mean_val', step_mean_val)
                    thisExp.addData('step_std_val', step_std_val)
                    thisExp.addData('block_number', block_number)
                    thisExp.addData('trial_number', trial_number)
                    thisExp.addData('correct', correct_ans)
                
                    # Increment trial number for next time
                    trial_number+=1
                    
                    # Check if we have any steps left
                    if steps_count >= steps_max:
                        steps_count = steps_max
                        predefined_setup_loop.finished = True
                        continueRoutine = False
                pts_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
                pts_polygon_abscissa.setSize((EDGE_DISTANCE*2,LINE_HEIGHT))
                pts_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
                pts_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
                pts_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
                pts_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
                pts_polygon_left_bound.setPos((-predefined_bound_x,PREDEFINED_BOUND_Y))
                pts_polygon_left_bound.setSize((PREDEFINED_BOUND_WIDTH,PREDEFINED_BOUND_HEIGHT))
                pts_polygon_left_bound.setLineColor(PREDEFINED_BOUND_COLOR)
                pts_polygon_right_bound.setPos((predefined_bound_x,PREDEFINED_BOUND_Y))
                pts_polygon_right_bound.setSize((PREDEFINED_BOUND_WIDTH,PREDEFINED_BOUND_HEIGHT))
                pts_polygon_stay_left.setPos((-predefined_bound_x_previous,PREDEFINED_BOUND_MARKER_Y))
                pts_polygon_stay_left.setSize((PREDEFINED_BOUND_MARKER_SZ,PREDEFINED_BOUND_MARKER_SZ))
                pts_polygon_stay_right.setPos((predefined_bound_x_previous,PREDEFINED_BOUND_MARKER_Y))
                pts_polygon_stay_right.setSize((PREDEFINED_BOUND_MARKER_SZ,PREDEFINED_BOUND_MARKER_SZ))
                pts_image_pigeon.setPos((pigeon_shown_x,ORIGIN_Y))
                pts_image_pigeon.setSize((pigeon_flip, 0.1))
                pts_polygon_coins_bar.setPos((COINS_BAR_X,coins_bar_center))
                pts_polygon_coins_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
                pts_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
                pts_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
                pts_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
                pts_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
                pts_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
                pts_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
                pts_text_coins.setPos((COINS_TEXT_X,COINS_TEXT_Y))
                pts_text_coins.setText(coins_count)
                pts_text_coins_properties.setText(coins_string)
                pts_text_steps.setText(steps_max - steps_count)
                pts_text_steps_properties.setText(steps_string)
                pts_text_instruction.setText(instruction_string)
                # store start times for predefined_task_setup
                predefined_task_setup.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
                predefined_task_setup.tStart = globalClock.getTime(format='float')
                predefined_task_setup.status = STARTED
                thisExp.addData('predefined_task_setup.started', predefined_task_setup.tStart)
                predefined_task_setup.maxDuration = None
                # keep track of which components have finished
                predefined_task_setupComponents = predefined_task_setup.components
                for thisComponent in predefined_task_setup.components:
                    thisComponent.tStart = None
                    thisComponent.tStop = None
                    thisComponent.tStartRefresh = None
                    thisComponent.tStopRefresh = None
                    if hasattr(thisComponent, 'status'):
                        thisComponent.status = NOT_STARTED
                # reset timers
                t = 0
                _timeToFirstFrame = win.getFutureFlipTime(clock="now")
                frameN = -1
                
                # --- Run Routine "predefined_task_setup" ---
                # if trial has changed, end Routine now
                if isinstance(predefined_setup_loop, data.TrialHandler2) and thisPredefined_setup_loop.thisN != predefined_setup_loop.thisTrial.thisN:
                    continueRoutine = False
                predefined_task_setup.forceEnded = routineForceEnded = not continueRoutine
                while continueRoutine:
                    # get current time
                    t = routineTimer.getTime()
                    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                    # update/draw components on each frame
                    
                    # *pts_kbd_update* updates
                    waitOnFlip = False
                    
                    # if pts_kbd_update is starting this frame...
                    if pts_kbd_update.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_kbd_update.frameNStart = frameN  # exact frame index
                        pts_kbd_update.tStart = t  # local t and not account for scr refresh
                        pts_kbd_update.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_kbd_update, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_kbd_update.status = STARTED
                        # keyboard checking is just starting
                        waitOnFlip = True
                        win.callOnFlip(pts_kbd_update.clock.reset)  # t=0 on next screen flip
                    if pts_kbd_update.status == STARTED and not waitOnFlip:
                        theseKeys = pts_kbd_update.getKeys(keyList=['up','down','left','right','space', 'c'], ignoreKeys=["escape"], waitRelease=False)
                        _pts_kbd_update_allKeys.extend(theseKeys)
                        if len(_pts_kbd_update_allKeys):
                            pts_kbd_update.keys = _pts_kbd_update_allKeys[-1].name  # just the last key pressed
                            pts_kbd_update.rt = _pts_kbd_update_allKeys[-1].rt
                            pts_kbd_update.duration = _pts_kbd_update_allKeys[-1].duration
                            # a response ends the routine
                            continueRoutine = False
                    
                    # *pts_polygon_abscissa* updates
                    
                    # if pts_polygon_abscissa is starting this frame...
                    if pts_polygon_abscissa.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_abscissa.frameNStart = frameN  # exact frame index
                        pts_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                        pts_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_abscissa.status = STARTED
                        pts_polygon_abscissa.setAutoDraw(True)
                    
                    # if pts_polygon_abscissa is active this frame...
                    if pts_polygon_abscissa.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_polygon_ordinate* updates
                    
                    # if pts_polygon_ordinate is starting this frame...
                    if pts_polygon_ordinate.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_ordinate.frameNStart = frameN  # exact frame index
                        pts_polygon_ordinate.tStart = t  # local t and not account for scr refresh
                        pts_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_ordinate.status = STARTED
                        pts_polygon_ordinate.setAutoDraw(True)
                    
                    # if pts_polygon_ordinate is active this frame...
                    if pts_polygon_ordinate.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_polygon_petey_midpoint* updates
                    
                    # if pts_polygon_petey_midpoint is starting this frame...
                    if pts_polygon_petey_midpoint.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
                        pts_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
                        pts_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_petey_midpoint.status = STARTED
                        pts_polygon_petey_midpoint.setAutoDraw(True)
                    
                    # if pts_polygon_petey_midpoint is active this frame...
                    if pts_polygon_petey_midpoint.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_polygon_left_bound* updates
                    
                    # if pts_polygon_left_bound is starting this frame...
                    if pts_polygon_left_bound.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_left_bound.frameNStart = frameN  # exact frame index
                        pts_polygon_left_bound.tStart = t  # local t and not account for scr refresh
                        pts_polygon_left_bound.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_left_bound, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_left_bound.status = STARTED
                        pts_polygon_left_bound.setAutoDraw(True)
                    
                    # if pts_polygon_left_bound is active this frame...
                    if pts_polygon_left_bound.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_polygon_right_bound* updates
                    
                    # if pts_polygon_right_bound is starting this frame...
                    if pts_polygon_right_bound.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_right_bound.frameNStart = frameN  # exact frame index
                        pts_polygon_right_bound.tStart = t  # local t and not account for scr refresh
                        pts_polygon_right_bound.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_right_bound, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_right_bound.status = STARTED
                        pts_polygon_right_bound.setAutoDraw(True)
                    
                    # if pts_polygon_right_bound is active this frame...
                    if pts_polygon_right_bound.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_polygon_stay_left* updates
                    
                    # if pts_polygon_stay_left is starting this frame...
                    if pts_polygon_stay_left.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_stay_left.frameNStart = frameN  # exact frame index
                        pts_polygon_stay_left.tStart = t  # local t and not account for scr refresh
                        pts_polygon_stay_left.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_stay_left, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_stay_left.status = STARTED
                        pts_polygon_stay_left.setAutoDraw(True)
                    
                    # if pts_polygon_stay_left is active this frame...
                    if pts_polygon_stay_left.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_polygon_stay_right* updates
                    
                    # if pts_polygon_stay_right is starting this frame...
                    if pts_polygon_stay_right.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_stay_right.frameNStart = frameN  # exact frame index
                        pts_polygon_stay_right.tStart = t  # local t and not account for scr refresh
                        pts_polygon_stay_right.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_stay_right, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_stay_right.status = STARTED
                        pts_polygon_stay_right.setAutoDraw(True)
                    
                    # if pts_polygon_stay_right is active this frame...
                    if pts_polygon_stay_right.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_image_pigeon* updates
                    
                    # if pts_image_pigeon is starting this frame...
                    if pts_image_pigeon.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_image_pigeon.frameNStart = frameN  # exact frame index
                        pts_image_pigeon.tStart = t  # local t and not account for scr refresh
                        pts_image_pigeon.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_image_pigeon, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_image_pigeon.status = STARTED
                        pts_image_pigeon.setAutoDraw(True)
                    
                    # if pts_image_pigeon is active this frame...
                    if pts_image_pigeon.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_image_right_seeds* updates
                    
                    # if pts_image_right_seeds is starting this frame...
                    if pts_image_right_seeds.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_image_right_seeds.frameNStart = frameN  # exact frame index
                        pts_image_right_seeds.tStart = t  # local t and not account for scr refresh
                        pts_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_image_right_seeds.status = STARTED
                        pts_image_right_seeds.setAutoDraw(True)
                    
                    # if pts_image_right_seeds is active this frame...
                    if pts_image_right_seeds.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_image_left_seeds* updates
                    
                    # if pts_image_left_seeds is starting this frame...
                    if pts_image_left_seeds.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_image_left_seeds.frameNStart = frameN  # exact frame index
                        pts_image_left_seeds.tStart = t  # local t and not account for scr refresh
                        pts_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_image_left_seeds.status = STARTED
                        pts_image_left_seeds.setAutoDraw(True)
                    
                    # if pts_image_left_seeds is active this frame...
                    if pts_image_left_seeds.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_polygon_coins_bar* updates
                    
                    # if pts_polygon_coins_bar is starting this frame...
                    if pts_polygon_coins_bar.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_coins_bar.frameNStart = frameN  # exact frame index
                        pts_polygon_coins_bar.tStart = t  # local t and not account for scr refresh
                        pts_polygon_coins_bar.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_coins_bar, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_coins_bar.status = STARTED
                        pts_polygon_coins_bar.setAutoDraw(True)
                    
                    # if pts_polygon_coins_bar is active this frame...
                    if pts_polygon_coins_bar.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_polygon_steps_bar* updates
                    
                    # if pts_polygon_steps_bar is starting this frame...
                    if pts_polygon_steps_bar.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_steps_bar.frameNStart = frameN  # exact frame index
                        pts_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
                        pts_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_steps_bar.status = STARTED
                        pts_polygon_steps_bar.setAutoDraw(True)
                    
                    # if pts_polygon_steps_bar is active this frame...
                    if pts_polygon_steps_bar.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_polygon_status_bar_apex* updates
                    
                    # if pts_polygon_status_bar_apex is starting this frame...
                    if pts_polygon_status_bar_apex.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
                        pts_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
                        pts_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_status_bar_apex.status = STARTED
                        pts_polygon_status_bar_apex.setAutoDraw(True)
                    
                    # if pts_polygon_status_bar_apex is active this frame...
                    if pts_polygon_status_bar_apex.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_polygon_status_bar_base* updates
                    
                    # if pts_polygon_status_bar_base is starting this frame...
                    if pts_polygon_status_bar_base.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_polygon_status_bar_base.frameNStart = frameN  # exact frame index
                        pts_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
                        pts_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_polygon_status_bar_base.status = STARTED
                        pts_polygon_status_bar_base.setAutoDraw(True)
                    
                    # if pts_polygon_status_bar_base is active this frame...
                    if pts_polygon_status_bar_base.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_text_coins* updates
                    
                    # if pts_text_coins is starting this frame...
                    if pts_text_coins.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_text_coins.frameNStart = frameN  # exact frame index
                        pts_text_coins.tStart = t  # local t and not account for scr refresh
                        pts_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_text_coins, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_text_coins.status = STARTED
                        pts_text_coins.setAutoDraw(True)
                    
                    # if pts_text_coins is active this frame...
                    if pts_text_coins.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_text_coins_label* updates
                    
                    # if pts_text_coins_label is starting this frame...
                    if pts_text_coins_label.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_text_coins_label.frameNStart = frameN  # exact frame index
                        pts_text_coins_label.tStart = t  # local t and not account for scr refresh
                        pts_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_text_coins_label, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_text_coins_label.status = STARTED
                        pts_text_coins_label.setAutoDraw(True)
                    
                    # if pts_text_coins_label is active this frame...
                    if pts_text_coins_label.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_text_coins_properties* updates
                    
                    # if pts_text_coins_properties is starting this frame...
                    if pts_text_coins_properties.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_text_coins_properties.frameNStart = frameN  # exact frame index
                        pts_text_coins_properties.tStart = t  # local t and not account for scr refresh
                        pts_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_text_coins_properties.status = STARTED
                        pts_text_coins_properties.setAutoDraw(True)
                    
                    # if pts_text_coins_properties is active this frame...
                    if pts_text_coins_properties.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_text_steps* updates
                    
                    # if pts_text_steps is starting this frame...
                    if pts_text_steps.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_text_steps.frameNStart = frameN  # exact frame index
                        pts_text_steps.tStart = t  # local t and not account for scr refresh
                        pts_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_text_steps, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_text_steps.status = STARTED
                        pts_text_steps.setAutoDraw(True)
                    
                    # if pts_text_steps is active this frame...
                    if pts_text_steps.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_text_steps_label* updates
                    
                    # if pts_text_steps_label is starting this frame...
                    if pts_text_steps_label.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_text_steps_label.frameNStart = frameN  # exact frame index
                        pts_text_steps_label.tStart = t  # local t and not account for scr refresh
                        pts_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_text_steps_label, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_text_steps_label.status = STARTED
                        pts_text_steps_label.setAutoDraw(True)
                    
                    # if pts_text_steps_label is active this frame...
                    if pts_text_steps_label.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_text_steps_properties* updates
                    
                    # if pts_text_steps_properties is starting this frame...
                    if pts_text_steps_properties.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_text_steps_properties.frameNStart = frameN  # exact frame index
                        pts_text_steps_properties.tStart = t  # local t and not account for scr refresh
                        pts_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_text_steps_properties.status = STARTED
                        pts_text_steps_properties.setAutoDraw(True)
                    
                    # if pts_text_steps_properties is active this frame...
                    if pts_text_steps_properties.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_text_instruction* updates
                    
                    # if pts_text_instruction is starting this frame...
                    if pts_text_instruction.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_text_instruction.frameNStart = frameN  # exact frame index
                        pts_text_instruction.tStart = t  # local t and not account for scr refresh
                        pts_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_text_instruction, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_text_instruction.status = STARTED
                        pts_text_instruction.setAutoDraw(True)
                    
                    # if pts_text_instruction is active this frame...
                    if pts_text_instruction.status == STARTED:
                        # update params
                        pass
                    
                    # *pts_text_hint* updates
                    
                    # if pts_text_hint is starting this frame...
                    if pts_text_hint.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_text_hint.frameNStart = frameN  # exact frame index
                        pts_text_hint.tStart = t  # local t and not account for scr refresh
                        pts_text_hint.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_text_hint, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_text_hint.status = STARTED
                        pts_text_hint.setAutoDraw(True)
                    
                    # if pts_text_hint is active this frame...
                    if pts_text_hint.status == STARTED:
                        # update params
                        pts_text_hint.setPos((ORIGIN_X,MESSAGE_Y), log=False)
                    
                    # *pts_text_hint_2* updates
                    
                    # if pts_text_hint_2 is starting this frame...
                    if pts_text_hint_2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        pts_text_hint_2.frameNStart = frameN  # exact frame index
                        pts_text_hint_2.tStart = t  # local t and not account for scr refresh
                        pts_text_hint_2.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(pts_text_hint_2, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        pts_text_hint_2.status = STARTED
                        pts_text_hint_2.setAutoDraw(True)
                    
                    # if pts_text_hint_2 is active this frame...
                    if pts_text_hint_2.status == STARTED:
                        # update params
                        pts_text_hint_2.setPos((ORIGIN_X,0.2), log=False)
                    
                    # check for quit (typically the Esc key)
                    if defaultKeyboard.getKeys(keyList=["escape"]):
                        thisExp.status = FINISHED
                    if thisExp.status == FINISHED or endExpNow:
                        endExperiment(thisExp, win=win)
                        return
                    # pause experiment here if requested
                    if thisExp.status == PAUSED:
                        pauseExperiment(
                            thisExp=thisExp, 
                            win=win, 
                            timers=[routineTimer], 
                            playbackComponents=[]
                        )
                        # skip the frame we paused on
                        continue
                    
                    # check if all components have finished
                    if not continueRoutine:  # a component has requested a forced-end of Routine
                        predefined_task_setup.forceEnded = routineForceEnded = True
                        break
                    continueRoutine = False  # will revert to True if at least one component still running
                    for thisComponent in predefined_task_setup.components:
                        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                            continueRoutine = True
                            break  # at least one component has not yet finished
                    
                    # refresh the screen
                    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                        win.flip()
                
                # --- Ending Routine "predefined_task_setup" ---
                for thisComponent in predefined_task_setup.components:
                    if hasattr(thisComponent, "setAutoDraw"):
                        thisComponent.setAutoDraw(False)
                # store stop times for predefined_task_setup
                predefined_task_setup.tStop = globalClock.getTime(format='float')
                predefined_task_setup.tStopRefresh = tThisFlipGlobal
                thisExp.addData('predefined_task_setup.stopped', predefined_task_setup.tStop)
                # Run 'End Routine' code from pts_code_update
                # Parse keypress
                if task_type == 'predefined_commit' and pts_kbd_update.keys == 'c' and (steps_count > min_steps_to_commit):
                    # Commit to a bound!
                    thisExp.addData('predefined_bound_final', predefined_bound_x)
                    thisExp.addData('predefined_bound', predefined_bound_x)
                    committed_to_bound = True
                    predefined_setup_loop.finished = True
                elif pts_kbd_update.keys == 'space':
                    # Save bound and continue
                    thisExp.addData('predefined_bound', predefined_bound_x)
                    predefined_setup_loop.finished = True
                elif pts_kbd_update.keys == 'left' and predefined_bound_x <= EDGE_DISTANCE-PREDEFINED_BOUND_SMALL_DELTA:
                    # Small-step larger
                    predefined_bound_x += PREDEFINED_BOUND_SMALL_DELTA
                elif pts_kbd_update.keys == 'right' and predefined_bound_x >= PREDEFINED_BOUND_SMALL_DELTA:
                    # Small-step smaller
                    predefined_bound_x -= PREDEFINED_BOUND_SMALL_DELTA
                elif pts_kbd_update.keys == 'up' and predefined_bound_x <= EDGE_DISTANCE-PREDEFINED_BOUND_LARGE_DELTA:
                    # Large-step larger
                    predefined_bound_x += PREDEFINED_BOUND_LARGE_DELTA
                elif pts_kbd_update.keys == 'down' and predefined_bound_x >= PREDEFINED_BOUND_LARGE_DELTA:
                    # Large-step smaller
                    predefined_bound_x -= PREDEFINED_BOUND_LARGE_DELTA
                
                pts_kbd_update.clearEvents()
                # the Routine "predefined_task_setup" was not non-slip safe, so reset the non-slip timer
                routineTimer.reset()
            # completed 1000.0 repeats of 'predefined_setup_loop'
            
            
            # set up handler to look after randomisation of conditions etc
            predefined_move_loop = data.TrialHandler2(
                name='predefined_move_loop',
                nReps=STEPS_MAX_PER_TRIAL, 
                method='random', 
                extraInfo=expInfo, 
                originPath=-1, 
                trialList=[None], 
                seed=None, 
            )
            thisExp.addLoop(predefined_move_loop)  # add the loop to the experiment
            thisPredefined_move_loop = predefined_move_loop.trialList[0]  # so we can initialise stimuli with some values
            # abbreviate parameter names if possible (e.g. rgb = thisPredefined_move_loop.rgb)
            if thisPredefined_move_loop != None:
                for paramName in thisPredefined_move_loop:
                    globals()[paramName] = thisPredefined_move_loop[paramName]
            
            for thisPredefined_move_loop in predefined_move_loop:
                currentLoop = predefined_move_loop
                thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
                # abbreviate parameter names if possible (e.g. rgb = thisPredefined_move_loop.rgb)
                if thisPredefined_move_loop != None:
                    for paramName in thisPredefined_move_loop:
                        globals()[paramName] = thisPredefined_move_loop[paramName]
                
                # --- Prepare to start Routine "predefined_task_run" ---
                # create an object to store info about Routine predefined_task_run
                predefined_task_run = data.Routine(
                    name='predefined_task_run',
                    components=[ptr_kbd_response, ptr_polygon_abscissa, ptr_polygon_ordinate, ptr_polygon_petey_midpoint, ptr_polygon_left_bound, ptr_polygon_right_bound, ptr_image_pigeon, ptr_image_right_seeds, ptr_image_left_seeds, ptr_polygon_coins_bar, ptr_polygon_steps_bar, ptr_polygon_status_bar_apex, ptr_polygon_status_bar_base, ptr_text_coins, ptr_text_coins_label, ptr_text_coins_properties, ptr_text_steps, ptr_text_steps_label, ptr_text_steps_properties, ptr_text_instruction],
                )
                predefined_task_run.status = NOT_STARTED
                continueRoutine = True
                # update component parameters for each repeat
                # create starting attributes for ptr_kbd_response
                ptr_kbd_response.keys = []
                ptr_kbd_response.rt = []
                _ptr_kbd_response_allKeys = []
                # Run 'Begin Routine' code from ptr_code_update
                # Check for end of task (no more steps)
                if (committed_to_bound and (steps_count > min_steps_to_commit)) or steps_count >= steps_max:
                    
                    # Clean up
                    if steps_count > steps_max:
                        steps_count = steps_max
                    predefined_move_loop.finished = True
                    predefined_trial_loop.finished = True
                    committed_to_bound = False
                    
                     # Show no steps remaining and end block string
                    steps_count = steps_max    
                    update_status(False,False,True)
                    set_block_end_string()
                    
                    # Save results    
                    thisExp.addData('steps_count', steps_count)    
                    thisExp.addData('pigeon_steps', pigeon_steps)
                    thisExp.addData('coins_count', coins_count)     
                    thisExp.addData('choice', 'none')
                    thisExp.addData('predefined_bound', predefined_bound_x)
                    thisExp.addData('predefined_bound_final', predefined_bound_x)
                    
                    # Set end of block message, which also determines bonus
                    set_block_end_string()
                    thisExp.addData('bonus', bonus_count)   
                
                else:
                    
                    # Increment the counter and move the pigeon
                    steps_count += 1
                    move_pigeon(step_mean_val, step_std_val)
                    pigeon_steps.append(pigeon_shown_x)
                
                    # Update steps status bar
                    update_status(False, False, True)
                
                    # Check bound crossings
                    still_going = False
                    if pigeon_true_x >= predefined_bound_x:
                        choice = 'right'
                    elif pigeon_true_x <= -predefined_bound_x:
                        choice = 'left'
                    elif steps_count >= steps_max:
                        choice = 'max'
                    else:        
                        still_going = True
                        
                    # Handle bound crossings
                    if not still_going:
                        
                        if choice == correct_ans:
                            coins_count += coins_gained_per_correct
                            set_trial_end_string(CORRECT_STRING)
                        elif choice != 'max':
                            coins_count -= coins_lost_per_error
                            steps_count += steps_lost_per_error
                            set_trial_end_string(ERROR_STRING)
                    
                        # Update coins status bar, feedback message
                        update_status(False, True, False)
                        
                        # Save results
                        thisExp.addData('pigeon_steps', pigeon_steps)
                        thisExp.addData('coins_count', coins_count)
                        thisExp.addData('steps_count', steps_count)
                        thisExp.addData('choice', choice)
                
                        # Drop out of pigeon movement loop
                        predefined_move_loop.finished = True
                
                ptr_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
                ptr_polygon_abscissa.setSize((EDGE_DISTANCE*2,LINE_HEIGHT))
                ptr_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
                ptr_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
                ptr_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
                ptr_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
                ptr_polygon_left_bound.setPos((-predefined_bound_x,PREDEFINED_BOUND_Y))
                ptr_polygon_left_bound.setSize((PREDEFINED_BOUND_WIDTH,PREDEFINED_BOUND_HEIGHT))
                ptr_polygon_left_bound.setLineColor(PREDEFINED_BOUND_COLOR)
                ptr_polygon_right_bound.setPos((predefined_bound_x,PREDEFINED_BOUND_Y))
                ptr_polygon_right_bound.setSize((PREDEFINED_BOUND_WIDTH,PREDEFINED_BOUND_HEIGHT))
                ptr_polygon_right_bound.setLineColor(PREDEFINED_BOUND_COLOR)
                ptr_image_pigeon.setPos((pigeon_shown_x,ORIGIN_Y))
                ptr_image_pigeon.setSize((pigeon_flip, 0.1))
                ptr_polygon_coins_bar.setPos((COINS_BAR_X,coins_bar_center))
                ptr_polygon_coins_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
                ptr_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
                ptr_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
                ptr_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
                ptr_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
                ptr_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
                ptr_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
                ptr_text_coins.setText(coins_count)
                ptr_text_coins_properties.setText(coins_string)
                ptr_text_steps.setText(steps_max - steps_count)
                ptr_text_steps_properties.setText(steps_string)
                ptr_text_instruction.setText(PREDEFINED_ABORT_INSTRUCTION)
                # store start times for predefined_task_run
                predefined_task_run.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
                predefined_task_run.tStart = globalClock.getTime(format='float')
                predefined_task_run.status = STARTED
                thisExp.addData('predefined_task_run.started', predefined_task_run.tStart)
                predefined_task_run.maxDuration = None
                # keep track of which components have finished
                predefined_task_runComponents = predefined_task_run.components
                for thisComponent in predefined_task_run.components:
                    thisComponent.tStart = None
                    thisComponent.tStop = None
                    thisComponent.tStartRefresh = None
                    thisComponent.tStopRefresh = None
                    if hasattr(thisComponent, 'status'):
                        thisComponent.status = NOT_STARTED
                # reset timers
                t = 0
                _timeToFirstFrame = win.getFutureFlipTime(clock="now")
                frameN = -1
                
                # --- Run Routine "predefined_task_run" ---
                # if trial has changed, end Routine now
                if isinstance(predefined_move_loop, data.TrialHandler2) and thisPredefined_move_loop.thisN != predefined_move_loop.thisTrial.thisN:
                    continueRoutine = False
                predefined_task_run.forceEnded = routineForceEnded = not continueRoutine
                while continueRoutine:
                    # get current time
                    t = routineTimer.getTime()
                    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                    # update/draw components on each frame
                    
                    # *ptr_kbd_response* updates
                    waitOnFlip = False
                    
                    # if ptr_kbd_response is starting this frame...
                    if ptr_kbd_response.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_kbd_response.frameNStart = frameN  # exact frame index
                        ptr_kbd_response.tStart = t  # local t and not account for scr refresh
                        ptr_kbd_response.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_kbd_response, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_kbd_response.status = STARTED
                        # keyboard checking is just starting
                        waitOnFlip = True
                        win.callOnFlip(ptr_kbd_response.clock.reset)  # t=0 on next screen flip
                    
                    # if ptr_kbd_response is stopping this frame...
                    if ptr_kbd_response.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_kbd_response.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_kbd_response.tStop = t  # not accounting for scr refresh
                            ptr_kbd_response.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_kbd_response.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_kbd_response.status = FINISHED
                            ptr_kbd_response.status = FINISHED
                    if ptr_kbd_response.status == STARTED and not waitOnFlip:
                        theseKeys = ptr_kbd_response.getKeys(keyList=['space', 'up', 'down', 'left', 'right'], ignoreKeys=["escape"], waitRelease=False)
                        _ptr_kbd_response_allKeys.extend(theseKeys)
                        if len(_ptr_kbd_response_allKeys):
                            ptr_kbd_response.keys = _ptr_kbd_response_allKeys[-1].name  # just the last key pressed
                            ptr_kbd_response.rt = _ptr_kbd_response_allKeys[-1].rt
                            ptr_kbd_response.duration = _ptr_kbd_response_allKeys[-1].duration
                            # a response ends the routine
                            continueRoutine = False
                    
                    # *ptr_polygon_abscissa* updates
                    
                    # if ptr_polygon_abscissa is starting this frame...
                    if ptr_polygon_abscissa.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_polygon_abscissa.frameNStart = frameN  # exact frame index
                        ptr_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                        ptr_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_polygon_abscissa.status = STARTED
                        ptr_polygon_abscissa.setAutoDraw(True)
                    
                    # if ptr_polygon_abscissa is active this frame...
                    if ptr_polygon_abscissa.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_polygon_abscissa is stopping this frame...
                    if ptr_polygon_abscissa.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_polygon_abscissa.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_polygon_abscissa.tStop = t  # not accounting for scr refresh
                            ptr_polygon_abscissa.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_polygon_abscissa.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_polygon_abscissa.status = FINISHED
                            ptr_polygon_abscissa.setAutoDraw(False)
                    
                    # *ptr_polygon_ordinate* updates
                    
                    # if ptr_polygon_ordinate is starting this frame...
                    if ptr_polygon_ordinate.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_polygon_ordinate.frameNStart = frameN  # exact frame index
                        ptr_polygon_ordinate.tStart = t  # local t and not account for scr refresh
                        ptr_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_polygon_ordinate.status = STARTED
                        ptr_polygon_ordinate.setAutoDraw(True)
                    
                    # if ptr_polygon_ordinate is active this frame...
                    if ptr_polygon_ordinate.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_polygon_ordinate is stopping this frame...
                    if ptr_polygon_ordinate.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_polygon_ordinate.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_polygon_ordinate.tStop = t  # not accounting for scr refresh
                            ptr_polygon_ordinate.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_polygon_ordinate.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_polygon_ordinate.status = FINISHED
                            ptr_polygon_ordinate.setAutoDraw(False)
                    
                    # *ptr_polygon_petey_midpoint* updates
                    
                    # if ptr_polygon_petey_midpoint is starting this frame...
                    if ptr_polygon_petey_midpoint.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
                        ptr_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
                        ptr_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_polygon_petey_midpoint.status = STARTED
                        ptr_polygon_petey_midpoint.setAutoDraw(True)
                    
                    # if ptr_polygon_petey_midpoint is active this frame...
                    if ptr_polygon_petey_midpoint.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_polygon_petey_midpoint is stopping this frame...
                    if ptr_polygon_petey_midpoint.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_polygon_petey_midpoint.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_polygon_petey_midpoint.tStop = t  # not accounting for scr refresh
                            ptr_polygon_petey_midpoint.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_polygon_petey_midpoint.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_polygon_petey_midpoint.status = FINISHED
                            ptr_polygon_petey_midpoint.setAutoDraw(False)
                    
                    # *ptr_polygon_left_bound* updates
                    
                    # if ptr_polygon_left_bound is starting this frame...
                    if ptr_polygon_left_bound.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_polygon_left_bound.frameNStart = frameN  # exact frame index
                        ptr_polygon_left_bound.tStart = t  # local t and not account for scr refresh
                        ptr_polygon_left_bound.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_polygon_left_bound, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_polygon_left_bound.status = STARTED
                        ptr_polygon_left_bound.setAutoDraw(True)
                    
                    # if ptr_polygon_left_bound is active this frame...
                    if ptr_polygon_left_bound.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_polygon_left_bound is stopping this frame...
                    if ptr_polygon_left_bound.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_polygon_left_bound.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_polygon_left_bound.tStop = t  # not accounting for scr refresh
                            ptr_polygon_left_bound.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_polygon_left_bound.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_polygon_left_bound.status = FINISHED
                            ptr_polygon_left_bound.setAutoDraw(False)
                    
                    # *ptr_polygon_right_bound* updates
                    
                    # if ptr_polygon_right_bound is starting this frame...
                    if ptr_polygon_right_bound.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_polygon_right_bound.frameNStart = frameN  # exact frame index
                        ptr_polygon_right_bound.tStart = t  # local t and not account for scr refresh
                        ptr_polygon_right_bound.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_polygon_right_bound, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_polygon_right_bound.status = STARTED
                        ptr_polygon_right_bound.setAutoDraw(True)
                    
                    # if ptr_polygon_right_bound is active this frame...
                    if ptr_polygon_right_bound.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_polygon_right_bound is stopping this frame...
                    if ptr_polygon_right_bound.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_polygon_right_bound.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_polygon_right_bound.tStop = t  # not accounting for scr refresh
                            ptr_polygon_right_bound.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_polygon_right_bound.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_polygon_right_bound.status = FINISHED
                            ptr_polygon_right_bound.setAutoDraw(False)
                    
                    # *ptr_image_pigeon* updates
                    
                    # if ptr_image_pigeon is starting this frame...
                    if ptr_image_pigeon.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_image_pigeon.frameNStart = frameN  # exact frame index
                        ptr_image_pigeon.tStart = t  # local t and not account for scr refresh
                        ptr_image_pigeon.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_image_pigeon, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_image_pigeon.status = STARTED
                        ptr_image_pigeon.setAutoDraw(True)
                    
                    # if ptr_image_pigeon is active this frame...
                    if ptr_image_pigeon.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_image_pigeon is stopping this frame...
                    if ptr_image_pigeon.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_image_pigeon.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_image_pigeon.tStop = t  # not accounting for scr refresh
                            ptr_image_pigeon.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_image_pigeon.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_image_pigeon.status = FINISHED
                            ptr_image_pigeon.setAutoDraw(False)
                    
                    # *ptr_image_right_seeds* updates
                    
                    # if ptr_image_right_seeds is starting this frame...
                    if ptr_image_right_seeds.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_image_right_seeds.frameNStart = frameN  # exact frame index
                        ptr_image_right_seeds.tStart = t  # local t and not account for scr refresh
                        ptr_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_image_right_seeds.status = STARTED
                        ptr_image_right_seeds.setAutoDraw(True)
                    
                    # if ptr_image_right_seeds is active this frame...
                    if ptr_image_right_seeds.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_image_right_seeds is stopping this frame...
                    if ptr_image_right_seeds.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_image_right_seeds.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_image_right_seeds.tStop = t  # not accounting for scr refresh
                            ptr_image_right_seeds.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_image_right_seeds.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_image_right_seeds.status = FINISHED
                            ptr_image_right_seeds.setAutoDraw(False)
                    
                    # *ptr_image_left_seeds* updates
                    
                    # if ptr_image_left_seeds is starting this frame...
                    if ptr_image_left_seeds.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_image_left_seeds.frameNStart = frameN  # exact frame index
                        ptr_image_left_seeds.tStart = t  # local t and not account for scr refresh
                        ptr_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_image_left_seeds.status = STARTED
                        ptr_image_left_seeds.setAutoDraw(True)
                    
                    # if ptr_image_left_seeds is active this frame...
                    if ptr_image_left_seeds.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_image_left_seeds is stopping this frame...
                    if ptr_image_left_seeds.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_image_left_seeds.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_image_left_seeds.tStop = t  # not accounting for scr refresh
                            ptr_image_left_seeds.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_image_left_seeds.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_image_left_seeds.status = FINISHED
                            ptr_image_left_seeds.setAutoDraw(False)
                    
                    # *ptr_polygon_coins_bar* updates
                    
                    # if ptr_polygon_coins_bar is starting this frame...
                    if ptr_polygon_coins_bar.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_polygon_coins_bar.frameNStart = frameN  # exact frame index
                        ptr_polygon_coins_bar.tStart = t  # local t and not account for scr refresh
                        ptr_polygon_coins_bar.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_polygon_coins_bar, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_polygon_coins_bar.status = STARTED
                        ptr_polygon_coins_bar.setAutoDraw(True)
                    
                    # if ptr_polygon_coins_bar is active this frame...
                    if ptr_polygon_coins_bar.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_polygon_coins_bar is stopping this frame...
                    if ptr_polygon_coins_bar.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_polygon_coins_bar.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_polygon_coins_bar.tStop = t  # not accounting for scr refresh
                            ptr_polygon_coins_bar.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_polygon_coins_bar.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_polygon_coins_bar.status = FINISHED
                            ptr_polygon_coins_bar.setAutoDraw(False)
                    
                    # *ptr_polygon_steps_bar* updates
                    
                    # if ptr_polygon_steps_bar is starting this frame...
                    if ptr_polygon_steps_bar.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_polygon_steps_bar.frameNStart = frameN  # exact frame index
                        ptr_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
                        ptr_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_polygon_steps_bar.status = STARTED
                        ptr_polygon_steps_bar.setAutoDraw(True)
                    
                    # if ptr_polygon_steps_bar is active this frame...
                    if ptr_polygon_steps_bar.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_polygon_steps_bar is stopping this frame...
                    if ptr_polygon_steps_bar.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_polygon_steps_bar.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_polygon_steps_bar.tStop = t  # not accounting for scr refresh
                            ptr_polygon_steps_bar.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_polygon_steps_bar.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_polygon_steps_bar.status = FINISHED
                            ptr_polygon_steps_bar.setAutoDraw(False)
                    
                    # *ptr_polygon_status_bar_apex* updates
                    
                    # if ptr_polygon_status_bar_apex is starting this frame...
                    if ptr_polygon_status_bar_apex.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
                        ptr_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
                        ptr_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_polygon_status_bar_apex.status = STARTED
                        ptr_polygon_status_bar_apex.setAutoDraw(True)
                    
                    # if ptr_polygon_status_bar_apex is active this frame...
                    if ptr_polygon_status_bar_apex.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_polygon_status_bar_apex is stopping this frame...
                    if ptr_polygon_status_bar_apex.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_polygon_status_bar_apex.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_polygon_status_bar_apex.tStop = t  # not accounting for scr refresh
                            ptr_polygon_status_bar_apex.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_polygon_status_bar_apex.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_polygon_status_bar_apex.status = FINISHED
                            ptr_polygon_status_bar_apex.setAutoDraw(False)
                    
                    # *ptr_polygon_status_bar_base* updates
                    
                    # if ptr_polygon_status_bar_base is starting this frame...
                    if ptr_polygon_status_bar_base.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_polygon_status_bar_base.frameNStart = frameN  # exact frame index
                        ptr_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
                        ptr_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_polygon_status_bar_base.status = STARTED
                        ptr_polygon_status_bar_base.setAutoDraw(True)
                    
                    # if ptr_polygon_status_bar_base is active this frame...
                    if ptr_polygon_status_bar_base.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_polygon_status_bar_base is stopping this frame...
                    if ptr_polygon_status_bar_base.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_polygon_status_bar_base.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_polygon_status_bar_base.tStop = t  # not accounting for scr refresh
                            ptr_polygon_status_bar_base.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_polygon_status_bar_base.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_polygon_status_bar_base.status = FINISHED
                            ptr_polygon_status_bar_base.setAutoDraw(False)
                    
                    # *ptr_text_coins* updates
                    
                    # if ptr_text_coins is starting this frame...
                    if ptr_text_coins.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_text_coins.frameNStart = frameN  # exact frame index
                        ptr_text_coins.tStart = t  # local t and not account for scr refresh
                        ptr_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_text_coins, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_text_coins.status = STARTED
                        ptr_text_coins.setAutoDraw(True)
                    
                    # if ptr_text_coins is active this frame...
                    if ptr_text_coins.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_text_coins is stopping this frame...
                    if ptr_text_coins.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_text_coins.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_text_coins.tStop = t  # not accounting for scr refresh
                            ptr_text_coins.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_text_coins.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_text_coins.status = FINISHED
                            ptr_text_coins.setAutoDraw(False)
                    
                    # *ptr_text_coins_label* updates
                    
                    # if ptr_text_coins_label is starting this frame...
                    if ptr_text_coins_label.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_text_coins_label.frameNStart = frameN  # exact frame index
                        ptr_text_coins_label.tStart = t  # local t and not account for scr refresh
                        ptr_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_text_coins_label, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_text_coins_label.status = STARTED
                        ptr_text_coins_label.setAutoDraw(True)
                    
                    # if ptr_text_coins_label is active this frame...
                    if ptr_text_coins_label.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_text_coins_label is stopping this frame...
                    if ptr_text_coins_label.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_text_coins_label.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_text_coins_label.tStop = t  # not accounting for scr refresh
                            ptr_text_coins_label.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_text_coins_label.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_text_coins_label.status = FINISHED
                            ptr_text_coins_label.setAutoDraw(False)
                    
                    # *ptr_text_coins_properties* updates
                    
                    # if ptr_text_coins_properties is starting this frame...
                    if ptr_text_coins_properties.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_text_coins_properties.frameNStart = frameN  # exact frame index
                        ptr_text_coins_properties.tStart = t  # local t and not account for scr refresh
                        ptr_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_text_coins_properties.status = STARTED
                        ptr_text_coins_properties.setAutoDraw(True)
                    
                    # if ptr_text_coins_properties is active this frame...
                    if ptr_text_coins_properties.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_text_coins_properties is stopping this frame...
                    if ptr_text_coins_properties.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_text_coins_properties.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_text_coins_properties.tStop = t  # not accounting for scr refresh
                            ptr_text_coins_properties.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_text_coins_properties.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_text_coins_properties.status = FINISHED
                            ptr_text_coins_properties.setAutoDraw(False)
                    
                    # *ptr_text_steps* updates
                    
                    # if ptr_text_steps is starting this frame...
                    if ptr_text_steps.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_text_steps.frameNStart = frameN  # exact frame index
                        ptr_text_steps.tStart = t  # local t and not account for scr refresh
                        ptr_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_text_steps, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_text_steps.status = STARTED
                        ptr_text_steps.setAutoDraw(True)
                    
                    # if ptr_text_steps is active this frame...
                    if ptr_text_steps.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_text_steps is stopping this frame...
                    if ptr_text_steps.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_text_steps.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_text_steps.tStop = t  # not accounting for scr refresh
                            ptr_text_steps.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_text_steps.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_text_steps.status = FINISHED
                            ptr_text_steps.setAutoDraw(False)
                    
                    # *ptr_text_steps_label* updates
                    
                    # if ptr_text_steps_label is starting this frame...
                    if ptr_text_steps_label.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_text_steps_label.frameNStart = frameN  # exact frame index
                        ptr_text_steps_label.tStart = t  # local t and not account for scr refresh
                        ptr_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_text_steps_label, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_text_steps_label.status = STARTED
                        ptr_text_steps_label.setAutoDraw(True)
                    
                    # if ptr_text_steps_label is active this frame...
                    if ptr_text_steps_label.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_text_steps_label is stopping this frame...
                    if ptr_text_steps_label.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_text_steps_label.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_text_steps_label.tStop = t  # not accounting for scr refresh
                            ptr_text_steps_label.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_text_steps_label.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_text_steps_label.status = FINISHED
                            ptr_text_steps_label.setAutoDraw(False)
                    
                    # *ptr_text_steps_properties* updates
                    
                    # if ptr_text_steps_properties is starting this frame...
                    if ptr_text_steps_properties.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_text_steps_properties.frameNStart = frameN  # exact frame index
                        ptr_text_steps_properties.tStart = t  # local t and not account for scr refresh
                        ptr_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_text_steps_properties.status = STARTED
                        ptr_text_steps_properties.setAutoDraw(True)
                    
                    # if ptr_text_steps_properties is active this frame...
                    if ptr_text_steps_properties.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_text_steps_properties is stopping this frame...
                    if ptr_text_steps_properties.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_text_steps_properties.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_text_steps_properties.tStop = t  # not accounting for scr refresh
                            ptr_text_steps_properties.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_text_steps_properties.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_text_steps_properties.status = FINISHED
                            ptr_text_steps_properties.setAutoDraw(False)
                    
                    # *ptr_text_instruction* updates
                    
                    # if ptr_text_instruction is starting this frame...
                    if ptr_text_instruction.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        ptr_text_instruction.frameNStart = frameN  # exact frame index
                        ptr_text_instruction.tStart = t  # local t and not account for scr refresh
                        ptr_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(ptr_text_instruction, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        ptr_text_instruction.status = STARTED
                        ptr_text_instruction.setAutoDraw(True)
                    
                    # if ptr_text_instruction is active this frame...
                    if ptr_text_instruction.status == STARTED:
                        # update params
                        pass
                    
                    # if ptr_text_instruction is stopping this frame...
                    if ptr_text_instruction.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > ptr_text_instruction.tStartRefresh + PREDEFINED_UPDATE_INTERVAL-frameTolerance:
                            # keep track of stop time/frame for later
                            ptr_text_instruction.tStop = t  # not accounting for scr refresh
                            ptr_text_instruction.tStopRefresh = tThisFlipGlobal  # on global time
                            ptr_text_instruction.frameNStop = frameN  # exact frame index
                            # update status
                            ptr_text_instruction.status = FINISHED
                            ptr_text_instruction.setAutoDraw(False)
                    
                    # check for quit (typically the Esc key)
                    if defaultKeyboard.getKeys(keyList=["escape"]):
                        thisExp.status = FINISHED
                    if thisExp.status == FINISHED or endExpNow:
                        endExperiment(thisExp, win=win)
                        return
                    # pause experiment here if requested
                    if thisExp.status == PAUSED:
                        pauseExperiment(
                            thisExp=thisExp, 
                            win=win, 
                            timers=[routineTimer], 
                            playbackComponents=[]
                        )
                        # skip the frame we paused on
                        continue
                    
                    # check if all components have finished
                    if not continueRoutine:  # a component has requested a forced-end of Routine
                        predefined_task_run.forceEnded = routineForceEnded = True
                        break
                    continueRoutine = False  # will revert to True if at least one component still running
                    for thisComponent in predefined_task_run.components:
                        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                            continueRoutine = True
                            break  # at least one component has not yet finished
                    
                    # refresh the screen
                    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                        win.flip()
                
                # --- Ending Routine "predefined_task_run" ---
                for thisComponent in predefined_task_run.components:
                    if hasattr(thisComponent, "setAutoDraw"):
                        thisComponent.setAutoDraw(False)
                # store stop times for predefined_task_run
                predefined_task_run.tStop = globalClock.getTime(format='float')
                predefined_task_run.tStopRefresh = tThisFlipGlobal
                thisExp.addData('predefined_task_run.stopped', predefined_task_run.tStop)
                # Run 'End Routine' code from ptr_code_update
                if ptr_kbd_response.keys:
                    
                    # Clear keyboard
                    ptr_kbd_response.clearEvents() 
                    
                    # Set no choice
                    choice = 'none'
                
                    # Set feedback message, update status, and abort
                    set_trial_end_string(PREDEFINED_ABORT_STRING)
                    predefined_move_loop.finished = True
                    
                # Save results
                thisExp.addData('pigeon_steps', pigeon_steps)
                thisExp.addData('coins_count', coins_count)
                thisExp.addData('steps_count', steps_count)
                thisExp.addData('choice', choice)
                # the Routine "predefined_task_run" was not non-slip safe, so reset the non-slip timer
                routineTimer.reset()
            # completed STEPS_MAX_PER_TRIAL repeats of 'predefined_move_loop'
            
            
            # --- Prepare to start Routine "show_message" ---
            # create an object to store info about Routine show_message
            show_message = data.Routine(
                name='show_message',
                components=[sm_kbd_continue, sm_polygon_box, sm_polygon_abscissa, sm_polygon_ordinate, sm_polygon_petey_midpoint, sm_image_right_seeds, sm_image_left_seeds, sm_image_petey, sm_polygon_coin_bar, sm_polygon_steps_bar, sm_polygon_status_bar_apex, sm_polygon_status_bar_base, sm_text_coins, sm_text_coins_label, sm_text_coins_properties, sm_text_steps, sm_text_steps_label, sm_text_steps_properties, sm_text_message, sm_text_instruction],
            )
            show_message.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            # create starting attributes for sm_kbd_continue
            sm_kbd_continue.keys = []
            sm_kbd_continue.rt = []
            _sm_kbd_continue_allKeys = []
            # Run 'Begin Routine' code from sm_code_update
            # Clear all keyboard events
            sm_kbd_continue.clearEvents()
            sm_polygon_box.setPos((box_x,box_y))
            sm_polygon_box.setSize((box_width,box_height))
            sm_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
            sm_polygon_abscissa.setSize((ABSCISSA_WIDTH,LINE_HEIGHT))
            sm_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
            sm_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
            sm_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
            sm_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
            sm_image_petey.setPos((pigeon_shown_x,ORIGIN_Y))
            sm_image_petey.setSize((pigeon_flip, 0.1))
            sm_polygon_coin_bar.setPos((COINS_BAR_X,coins_bar_center))
            sm_polygon_coin_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
            sm_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
            sm_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
            sm_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
            sm_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
            sm_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
            sm_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
            sm_text_coins.setText(coins_count)
            sm_text_coins_properties.setText(coins_string)
            sm_text_steps.setText(steps_max - steps_count)
            sm_text_steps_properties.setText(steps_string)
            sm_text_message.setText(message_string)
            sm_text_instruction.setText(CONTINUE_STRING)
            # store start times for show_message
            show_message.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            show_message.tStart = globalClock.getTime(format='float')
            show_message.status = STARTED
            thisExp.addData('show_message.started', show_message.tStart)
            show_message.maxDuration = None
            # keep track of which components have finished
            show_messageComponents = show_message.components
            for thisComponent in show_message.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "show_message" ---
            # if trial has changed, end Routine now
            if isinstance(predefined_trial_loop, data.TrialHandler2) and thisPredefined_trial_loop.thisN != predefined_trial_loop.thisTrial.thisN:
                continueRoutine = False
            show_message.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                # *sm_kbd_continue* updates
                waitOnFlip = False
                
                # if sm_kbd_continue is starting this frame...
                if sm_kbd_continue.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    sm_kbd_continue.frameNStart = frameN  # exact frame index
                    sm_kbd_continue.tStart = t  # local t and not account for scr refresh
                    sm_kbd_continue.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_kbd_continue, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_kbd_continue.status = STARTED
                    # keyboard checking is just starting
                    waitOnFlip = True
                    win.callOnFlip(sm_kbd_continue.clock.reset)  # t=0 on next screen flip
                    win.callOnFlip(sm_kbd_continue.clearEvents, eventType='keyboard')  # clear events on next screen flip
                if sm_kbd_continue.status == STARTED and not waitOnFlip:
                    theseKeys = sm_kbd_continue.getKeys(keyList=['space', 'up', 'down', 'left', 'right'], ignoreKeys=["escape"], waitRelease=False)
                    _sm_kbd_continue_allKeys.extend(theseKeys)
                    if len(_sm_kbd_continue_allKeys):
                        sm_kbd_continue.keys = _sm_kbd_continue_allKeys[-1].name  # just the last key pressed
                        sm_kbd_continue.rt = _sm_kbd_continue_allKeys[-1].rt
                        sm_kbd_continue.duration = _sm_kbd_continue_allKeys[-1].duration
                        # a response ends the routine
                        continueRoutine = False
                
                # *sm_polygon_box* updates
                
                # if sm_polygon_box is starting this frame...
                if sm_polygon_box.status == NOT_STARTED and box_show:
                    # keep track of start time/frame for later
                    sm_polygon_box.frameNStart = frameN  # exact frame index
                    sm_polygon_box.tStart = t  # local t and not account for scr refresh
                    sm_polygon_box.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_box, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_box.status = STARTED
                    sm_polygon_box.setAutoDraw(True)
                
                # if sm_polygon_box is active this frame...
                if sm_polygon_box.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_abscissa* updates
                
                # if sm_polygon_abscissa is starting this frame...
                if sm_polygon_abscissa.status == NOT_STARTED and pigeon_show:
                    # keep track of start time/frame for later
                    sm_polygon_abscissa.frameNStart = frameN  # exact frame index
                    sm_polygon_abscissa.tStart = t  # local t and not account for scr refresh
                    sm_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_abscissa.status = STARTED
                    sm_polygon_abscissa.setAutoDraw(True)
                
                # if sm_polygon_abscissa is active this frame...
                if sm_polygon_abscissa.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_ordinate* updates
                
                # if sm_polygon_ordinate is starting this frame...
                if sm_polygon_ordinate.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_ordinate.frameNStart = frameN  # exact frame index
                    sm_polygon_ordinate.tStart = t  # local t and not account for scr refresh
                    sm_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_ordinate.status = STARTED
                    sm_polygon_ordinate.setAutoDraw(True)
                
                # if sm_polygon_ordinate is active this frame...
                if sm_polygon_ordinate.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_petey_midpoint* updates
                
                # if sm_polygon_petey_midpoint is starting this frame...
                if sm_polygon_petey_midpoint.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
                    sm_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
                    sm_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_petey_midpoint.status = STARTED
                    sm_polygon_petey_midpoint.setAutoDraw(True)
                
                # if sm_polygon_petey_midpoint is active this frame...
                if sm_polygon_petey_midpoint.status == STARTED:
                    # update params
                    pass
                
                # *sm_image_right_seeds* updates
                
                # if sm_image_right_seeds is starting this frame...
                if sm_image_right_seeds.status == NOT_STARTED and seeds_show:
                    # keep track of start time/frame for later
                    sm_image_right_seeds.frameNStart = frameN  # exact frame index
                    sm_image_right_seeds.tStart = t  # local t and not account for scr refresh
                    sm_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_image_right_seeds.status = STARTED
                    sm_image_right_seeds.setAutoDraw(True)
                
                # if sm_image_right_seeds is active this frame...
                if sm_image_right_seeds.status == STARTED:
                    # update params
                    pass
                
                # *sm_image_left_seeds* updates
                
                # if sm_image_left_seeds is starting this frame...
                if sm_image_left_seeds.status == NOT_STARTED and seeds_show:
                    # keep track of start time/frame for later
                    sm_image_left_seeds.frameNStart = frameN  # exact frame index
                    sm_image_left_seeds.tStart = t  # local t and not account for scr refresh
                    sm_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_image_left_seeds.status = STARTED
                    sm_image_left_seeds.setAutoDraw(True)
                
                # if sm_image_left_seeds is active this frame...
                if sm_image_left_seeds.status == STARTED:
                    # update params
                    pass
                
                # *sm_image_petey* updates
                
                # if sm_image_petey is starting this frame...
                if sm_image_petey.status == NOT_STARTED and pigeon_show:
                    # keep track of start time/frame for later
                    sm_image_petey.frameNStart = frameN  # exact frame index
                    sm_image_petey.tStart = t  # local t and not account for scr refresh
                    sm_image_petey.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_image_petey, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_image_petey.status = STARTED
                    sm_image_petey.setAutoDraw(True)
                
                # if sm_image_petey is active this frame...
                if sm_image_petey.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_coin_bar* updates
                
                # if sm_polygon_coin_bar is starting this frame...
                if sm_polygon_coin_bar.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_coin_bar.frameNStart = frameN  # exact frame index
                    sm_polygon_coin_bar.tStart = t  # local t and not account for scr refresh
                    sm_polygon_coin_bar.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_coin_bar, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_coin_bar.status = STARTED
                    sm_polygon_coin_bar.setAutoDraw(True)
                
                # if sm_polygon_coin_bar is active this frame...
                if sm_polygon_coin_bar.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_steps_bar* updates
                
                # if sm_polygon_steps_bar is starting this frame...
                if sm_polygon_steps_bar.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_steps_bar.frameNStart = frameN  # exact frame index
                    sm_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
                    sm_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_steps_bar.status = STARTED
                    sm_polygon_steps_bar.setAutoDraw(True)
                
                # if sm_polygon_steps_bar is active this frame...
                if sm_polygon_steps_bar.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_status_bar_apex* updates
                
                # if sm_polygon_status_bar_apex is starting this frame...
                if sm_polygon_status_bar_apex.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
                    sm_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
                    sm_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_status_bar_apex.status = STARTED
                    sm_polygon_status_bar_apex.setAutoDraw(True)
                
                # if sm_polygon_status_bar_apex is active this frame...
                if sm_polygon_status_bar_apex.status == STARTED:
                    # update params
                    pass
                
                # *sm_polygon_status_bar_base* updates
                
                # if sm_polygon_status_bar_base is starting this frame...
                if sm_polygon_status_bar_base.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_polygon_status_bar_base.frameNStart = frameN  # exact frame index
                    sm_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
                    sm_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_polygon_status_bar_base.status = STARTED
                    sm_polygon_status_bar_base.setAutoDraw(True)
                
                # if sm_polygon_status_bar_base is active this frame...
                if sm_polygon_status_bar_base.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_coins* updates
                
                # if sm_text_coins is starting this frame...
                if sm_text_coins.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_coins.frameNStart = frameN  # exact frame index
                    sm_text_coins.tStart = t  # local t and not account for scr refresh
                    sm_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_coins, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_coins.status = STARTED
                    sm_text_coins.setAutoDraw(True)
                
                # if sm_text_coins is active this frame...
                if sm_text_coins.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_coins_label* updates
                
                # if sm_text_coins_label is starting this frame...
                if sm_text_coins_label.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_coins_label.frameNStart = frameN  # exact frame index
                    sm_text_coins_label.tStart = t  # local t and not account for scr refresh
                    sm_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_coins_label, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_coins_label.status = STARTED
                    sm_text_coins_label.setAutoDraw(True)
                
                # if sm_text_coins_label is active this frame...
                if sm_text_coins_label.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_coins_properties* updates
                
                # if sm_text_coins_properties is starting this frame...
                if sm_text_coins_properties.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_coins_properties.frameNStart = frameN  # exact frame index
                    sm_text_coins_properties.tStart = t  # local t and not account for scr refresh
                    sm_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_coins_properties.status = STARTED
                    sm_text_coins_properties.setAutoDraw(True)
                
                # if sm_text_coins_properties is active this frame...
                if sm_text_coins_properties.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_steps* updates
                
                # if sm_text_steps is starting this frame...
                if sm_text_steps.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_steps.frameNStart = frameN  # exact frame index
                    sm_text_steps.tStart = t  # local t and not account for scr refresh
                    sm_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_steps, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_steps.status = STARTED
                    sm_text_steps.setAutoDraw(True)
                
                # if sm_text_steps is active this frame...
                if sm_text_steps.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_steps_label* updates
                
                # if sm_text_steps_label is starting this frame...
                if sm_text_steps_label.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_steps_label.frameNStart = frameN  # exact frame index
                    sm_text_steps_label.tStart = t  # local t and not account for scr refresh
                    sm_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_steps_label, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_steps_label.status = STARTED
                    sm_text_steps_label.setAutoDraw(True)
                
                # if sm_text_steps_label is active this frame...
                if sm_text_steps_label.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_steps_properties* updates
                
                # if sm_text_steps_properties is starting this frame...
                if sm_text_steps_properties.status == NOT_STARTED and status_show:
                    # keep track of start time/frame for later
                    sm_text_steps_properties.frameNStart = frameN  # exact frame index
                    sm_text_steps_properties.tStart = t  # local t and not account for scr refresh
                    sm_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_steps_properties.status = STARTED
                    sm_text_steps_properties.setAutoDraw(True)
                
                # if sm_text_steps_properties is active this frame...
                if sm_text_steps_properties.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_message* updates
                
                # if sm_text_message is starting this frame...
                if sm_text_message.status == NOT_STARTED and message_show:
                    # keep track of start time/frame for later
                    sm_text_message.frameNStart = frameN  # exact frame index
                    sm_text_message.tStart = t  # local t and not account for scr refresh
                    sm_text_message.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_message, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_message.status = STARTED
                    sm_text_message.setAutoDraw(True)
                
                # if sm_text_message is active this frame...
                if sm_text_message.status == STARTED:
                    # update params
                    pass
                
                # *sm_text_instruction* updates
                
                # if sm_text_instruction is starting this frame...
                if sm_text_instruction.status == NOT_STARTED and instruction_show:
                    # keep track of start time/frame for later
                    sm_text_instruction.frameNStart = frameN  # exact frame index
                    sm_text_instruction.tStart = t  # local t and not account for scr refresh
                    sm_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(sm_text_instruction, 'tStartRefresh')  # time at next scr refresh
                    # update status
                    sm_text_instruction.status = STARTED
                    sm_text_instruction.setAutoDraw(True)
                
                # if sm_text_instruction is active this frame...
                if sm_text_instruction.status == STARTED:
                    # update params
                    pass
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer], 
                        playbackComponents=[]
                    )
                    # skip the frame we paused on
                    continue
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    show_message.forceEnded = routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in show_message.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "show_message" ---
            for thisComponent in show_message.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for show_message
            show_message.tStop = globalClock.getTime(format='float')
            show_message.tStopRefresh = tThisFlipGlobal
            thisExp.addData('show_message.stopped', show_message.tStop)
            # Run 'End Routine' code from sm_code_update
            # Ick, but whatever. Makes sure that we don't ever
            # drop unexpectedly out of a move loop without 
            # resetting
            trial_start = True
            # the Routine "show_message" was not non-slip safe, so reset the non-slip timer
            routineTimer.reset()
            thisExp.nextEntry()
            
        # completed predefined_trials_count repeats of 'predefined_trial_loop'
        
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        thisExp.nextEntry()
        
    # completed 1.0 repeats of 'experiment_loop'
    
    if thisSession is not None:
        # if running in a Session with a Liaison client, send data up to now
        thisSession.sendExperimentData()
    
    # --- Prepare to start Routine "final_message" ---
    # create an object to store info about Routine final_message
    final_message = data.Routine(
        name='final_message',
        components=[],
    )
    final_message.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # Run 'Begin Routine' code from fm_code_setup
    message_string = 'Finished, fantastic job!\nYou earned ' + str(bonus_count) + ' bonuses.'
    # store start times for final_message
    final_message.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    final_message.tStart = globalClock.getTime(format='float')
    final_message.status = STARTED
    thisExp.addData('final_message.started', final_message.tStart)
    final_message.maxDuration = None
    # keep track of which components have finished
    final_messageComponents = final_message.components
    for thisComponent in final_message.components:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "final_message" ---
    final_message.forceEnded = routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, win=win)
            return
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
            )
            # skip the frame we paused on
            continue
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            final_message.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in final_message.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "final_message" ---
    for thisComponent in final_message.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for final_message
    final_message.tStop = globalClock.getTime(format='float')
    final_message.tStopRefresh = tThisFlipGlobal
    thisExp.addData('final_message.stopped', final_message.tStop)
    thisExp.nextEntry()
    # the Routine "final_message" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "show_final_message" ---
    # create an object to store info about Routine show_final_message
    show_final_message = data.Routine(
        name='show_final_message',
        components=[sfm_kbd_continue, sfm_text_message, sfm_text_steps_properties, sfm_text_steps_label, sfm_text_steps, sfm_text_coins_properties, sfm_text_coins_label, sfm_text_coins, sfm_polygon_status_bar_base, sfm_polygon_status_bar_apex, sfm_polygon_steps_bar, sfm_polygon_coin_bar, sfm_image_petey, sfm_image_left_seeds, sfm_image_right_seeds, sfm_polygon_petey_midpoint, sfm_polygon_ordinate, sfm_polygon_abscissa, sfm_polygon_box, sfm_text_instruction],
    )
    show_final_message.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    # create starting attributes for sfm_kbd_continue
    sfm_kbd_continue.keys = []
    sfm_kbd_continue.rt = []
    _sfm_kbd_continue_allKeys = []
    # Run 'Begin Routine' code from sfm_code_update
    # Clear all keyboard events
    sm_kbd_continue.clearEvents()
    sfm_text_message.setText(message_string)
    sfm_text_steps_properties.setText(steps_string)
    sfm_text_steps.setText(steps_max - steps_count)
    sfm_text_coins_properties.setText(coins_string)
    sfm_text_coins.setText(coins_count)
    sfm_polygon_status_bar_base.setPos((ORIGIN_X,STATUS_BAR_BASE_Y))
    sfm_polygon_status_bar_base.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
    sfm_polygon_status_bar_apex.setPos((ORIGIN_X,STATUS_BAR_APEX_Y))
    sfm_polygon_status_bar_apex.setSize((STATUS_BAR_BRACKET_WIDTH,LINE_HEIGHT))
    sfm_polygon_steps_bar.setPos((STEPS_BAR_X,steps_bar_center))
    sfm_polygon_steps_bar.setSize((STATUS_BAR_WIDTH,steps_bar_height))
    sfm_polygon_coin_bar.setPos((COINS_BAR_X,coins_bar_center))
    sfm_polygon_coin_bar.setSize((STATUS_BAR_WIDTH,coins_bar_height))
    sfm_image_petey.setPos((pigeon_shown_x,ORIGIN_Y))
    sfm_image_petey.setSize((pigeon_flip, 0.1))
    sfm_polygon_petey_midpoint.setPos((pigeon_shown_x,ORDINATE_Y))
    sfm_polygon_petey_midpoint.setSize((ORDINATE_WIDTH*0.9,ORDINATE_HEIGHT))
    sfm_polygon_ordinate.setPos((ORIGIN_X,ORDINATE_Y))
    sfm_polygon_ordinate.setSize((ORDINATE_WIDTH,ORDINATE_HEIGHT))
    sfm_polygon_abscissa.setPos((ORIGIN_X,ABSCISSA_Y))
    sfm_polygon_abscissa.setSize((ABSCISSA_WIDTH,LINE_HEIGHT))
    sfm_polygon_box.setPos((box_x,box_y))
    sfm_polygon_box.setSize((box_width,box_height))
    sfm_text_instruction.setText(CONTINUE_STRING)
    # store start times for show_final_message
    show_final_message.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    show_final_message.tStart = globalClock.getTime(format='float')
    show_final_message.status = STARTED
    thisExp.addData('show_final_message.started', show_final_message.tStart)
    show_final_message.maxDuration = None
    # keep track of which components have finished
    show_final_messageComponents = show_final_message.components
    for thisComponent in show_final_message.components:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "show_final_message" ---
    show_final_message.forceEnded = routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *sfm_kbd_continue* updates
        waitOnFlip = False
        
        # if sfm_kbd_continue is starting this frame...
        if sfm_kbd_continue.status == NOT_STARTED and tThisFlip >= 0.5-frameTolerance:
            # keep track of start time/frame for later
            sfm_kbd_continue.frameNStart = frameN  # exact frame index
            sfm_kbd_continue.tStart = t  # local t and not account for scr refresh
            sfm_kbd_continue.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_kbd_continue, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_kbd_continue.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(sfm_kbd_continue.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(sfm_kbd_continue.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if sfm_kbd_continue.status == STARTED and not waitOnFlip:
            theseKeys = sfm_kbd_continue.getKeys(keyList=['space', 'up', 'down', 'left', 'right'], ignoreKeys=["escape"], waitRelease=False)
            _sfm_kbd_continue_allKeys.extend(theseKeys)
            if len(_sfm_kbd_continue_allKeys):
                sfm_kbd_continue.keys = _sfm_kbd_continue_allKeys[-1].name  # just the last key pressed
                sfm_kbd_continue.rt = _sfm_kbd_continue_allKeys[-1].rt
                sfm_kbd_continue.duration = _sfm_kbd_continue_allKeys[-1].duration
                # a response ends the routine
                continueRoutine = False
        
        # *sfm_text_message* updates
        
        # if sfm_text_message is starting this frame...
        if sfm_text_message.status == NOT_STARTED and message_show:
            # keep track of start time/frame for later
            sfm_text_message.frameNStart = frameN  # exact frame index
            sfm_text_message.tStart = t  # local t and not account for scr refresh
            sfm_text_message.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_text_message, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_text_message.status = STARTED
            sfm_text_message.setAutoDraw(True)
        
        # if sfm_text_message is active this frame...
        if sfm_text_message.status == STARTED:
            # update params
            pass
        
        # *sfm_text_steps_properties* updates
        
        # if sfm_text_steps_properties is starting this frame...
        if sfm_text_steps_properties.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_text_steps_properties.frameNStart = frameN  # exact frame index
            sfm_text_steps_properties.tStart = t  # local t and not account for scr refresh
            sfm_text_steps_properties.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_text_steps_properties, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_text_steps_properties.status = STARTED
            sfm_text_steps_properties.setAutoDraw(True)
        
        # if sfm_text_steps_properties is active this frame...
        if sfm_text_steps_properties.status == STARTED:
            # update params
            pass
        
        # *sfm_text_steps_label* updates
        
        # if sfm_text_steps_label is starting this frame...
        if sfm_text_steps_label.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_text_steps_label.frameNStart = frameN  # exact frame index
            sfm_text_steps_label.tStart = t  # local t and not account for scr refresh
            sfm_text_steps_label.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_text_steps_label, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_text_steps_label.status = STARTED
            sfm_text_steps_label.setAutoDraw(True)
        
        # if sfm_text_steps_label is active this frame...
        if sfm_text_steps_label.status == STARTED:
            # update params
            pass
        
        # *sfm_text_steps* updates
        
        # if sfm_text_steps is starting this frame...
        if sfm_text_steps.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_text_steps.frameNStart = frameN  # exact frame index
            sfm_text_steps.tStart = t  # local t and not account for scr refresh
            sfm_text_steps.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_text_steps, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_text_steps.status = STARTED
            sfm_text_steps.setAutoDraw(True)
        
        # if sfm_text_steps is active this frame...
        if sfm_text_steps.status == STARTED:
            # update params
            pass
        
        # *sfm_text_coins_properties* updates
        
        # if sfm_text_coins_properties is starting this frame...
        if sfm_text_coins_properties.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_text_coins_properties.frameNStart = frameN  # exact frame index
            sfm_text_coins_properties.tStart = t  # local t and not account for scr refresh
            sfm_text_coins_properties.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_text_coins_properties, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_text_coins_properties.status = STARTED
            sfm_text_coins_properties.setAutoDraw(True)
        
        # if sfm_text_coins_properties is active this frame...
        if sfm_text_coins_properties.status == STARTED:
            # update params
            pass
        
        # *sfm_text_coins_label* updates
        
        # if sfm_text_coins_label is starting this frame...
        if sfm_text_coins_label.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_text_coins_label.frameNStart = frameN  # exact frame index
            sfm_text_coins_label.tStart = t  # local t and not account for scr refresh
            sfm_text_coins_label.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_text_coins_label, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_text_coins_label.status = STARTED
            sfm_text_coins_label.setAutoDraw(True)
        
        # if sfm_text_coins_label is active this frame...
        if sfm_text_coins_label.status == STARTED:
            # update params
            pass
        
        # *sfm_text_coins* updates
        
        # if sfm_text_coins is starting this frame...
        if sfm_text_coins.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_text_coins.frameNStart = frameN  # exact frame index
            sfm_text_coins.tStart = t  # local t and not account for scr refresh
            sfm_text_coins.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_text_coins, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_text_coins.status = STARTED
            sfm_text_coins.setAutoDraw(True)
        
        # if sfm_text_coins is active this frame...
        if sfm_text_coins.status == STARTED:
            # update params
            pass
        
        # *sfm_polygon_status_bar_base* updates
        
        # if sfm_polygon_status_bar_base is starting this frame...
        if sfm_polygon_status_bar_base.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_polygon_status_bar_base.frameNStart = frameN  # exact frame index
            sfm_polygon_status_bar_base.tStart = t  # local t and not account for scr refresh
            sfm_polygon_status_bar_base.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_polygon_status_bar_base, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_polygon_status_bar_base.status = STARTED
            sfm_polygon_status_bar_base.setAutoDraw(True)
        
        # if sfm_polygon_status_bar_base is active this frame...
        if sfm_polygon_status_bar_base.status == STARTED:
            # update params
            pass
        
        # *sfm_polygon_status_bar_apex* updates
        
        # if sfm_polygon_status_bar_apex is starting this frame...
        if sfm_polygon_status_bar_apex.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_polygon_status_bar_apex.frameNStart = frameN  # exact frame index
            sfm_polygon_status_bar_apex.tStart = t  # local t and not account for scr refresh
            sfm_polygon_status_bar_apex.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_polygon_status_bar_apex, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_polygon_status_bar_apex.status = STARTED
            sfm_polygon_status_bar_apex.setAutoDraw(True)
        
        # if sfm_polygon_status_bar_apex is active this frame...
        if sfm_polygon_status_bar_apex.status == STARTED:
            # update params
            pass
        
        # *sfm_polygon_steps_bar* updates
        
        # if sfm_polygon_steps_bar is starting this frame...
        if sfm_polygon_steps_bar.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_polygon_steps_bar.frameNStart = frameN  # exact frame index
            sfm_polygon_steps_bar.tStart = t  # local t and not account for scr refresh
            sfm_polygon_steps_bar.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_polygon_steps_bar, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_polygon_steps_bar.status = STARTED
            sfm_polygon_steps_bar.setAutoDraw(True)
        
        # if sfm_polygon_steps_bar is active this frame...
        if sfm_polygon_steps_bar.status == STARTED:
            # update params
            pass
        
        # *sfm_polygon_coin_bar* updates
        
        # if sfm_polygon_coin_bar is starting this frame...
        if sfm_polygon_coin_bar.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_polygon_coin_bar.frameNStart = frameN  # exact frame index
            sfm_polygon_coin_bar.tStart = t  # local t and not account for scr refresh
            sfm_polygon_coin_bar.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_polygon_coin_bar, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_polygon_coin_bar.status = STARTED
            sfm_polygon_coin_bar.setAutoDraw(True)
        
        # if sfm_polygon_coin_bar is active this frame...
        if sfm_polygon_coin_bar.status == STARTED:
            # update params
            pass
        
        # *sfm_image_petey* updates
        
        # if sfm_image_petey is starting this frame...
        if sfm_image_petey.status == NOT_STARTED and pigeon_show:
            # keep track of start time/frame for later
            sfm_image_petey.frameNStart = frameN  # exact frame index
            sfm_image_petey.tStart = t  # local t and not account for scr refresh
            sfm_image_petey.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_image_petey, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_image_petey.status = STARTED
            sfm_image_petey.setAutoDraw(True)
        
        # if sfm_image_petey is active this frame...
        if sfm_image_petey.status == STARTED:
            # update params
            pass
        
        # *sfm_image_left_seeds* updates
        
        # if sfm_image_left_seeds is starting this frame...
        if sfm_image_left_seeds.status == NOT_STARTED and seeds_show:
            # keep track of start time/frame for later
            sfm_image_left_seeds.frameNStart = frameN  # exact frame index
            sfm_image_left_seeds.tStart = t  # local t and not account for scr refresh
            sfm_image_left_seeds.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_image_left_seeds, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_image_left_seeds.status = STARTED
            sfm_image_left_seeds.setAutoDraw(True)
        
        # if sfm_image_left_seeds is active this frame...
        if sfm_image_left_seeds.status == STARTED:
            # update params
            pass
        
        # *sfm_image_right_seeds* updates
        
        # if sfm_image_right_seeds is starting this frame...
        if sfm_image_right_seeds.status == NOT_STARTED and seeds_show:
            # keep track of start time/frame for later
            sfm_image_right_seeds.frameNStart = frameN  # exact frame index
            sfm_image_right_seeds.tStart = t  # local t and not account for scr refresh
            sfm_image_right_seeds.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_image_right_seeds, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_image_right_seeds.status = STARTED
            sfm_image_right_seeds.setAutoDraw(True)
        
        # if sfm_image_right_seeds is active this frame...
        if sfm_image_right_seeds.status == STARTED:
            # update params
            pass
        
        # *sfm_polygon_petey_midpoint* updates
        
        # if sfm_polygon_petey_midpoint is starting this frame...
        if sfm_polygon_petey_midpoint.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_polygon_petey_midpoint.frameNStart = frameN  # exact frame index
            sfm_polygon_petey_midpoint.tStart = t  # local t and not account for scr refresh
            sfm_polygon_petey_midpoint.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_polygon_petey_midpoint, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_polygon_petey_midpoint.status = STARTED
            sfm_polygon_petey_midpoint.setAutoDraw(True)
        
        # if sfm_polygon_petey_midpoint is active this frame...
        if sfm_polygon_petey_midpoint.status == STARTED:
            # update params
            pass
        
        # *sfm_polygon_ordinate* updates
        
        # if sfm_polygon_ordinate is starting this frame...
        if sfm_polygon_ordinate.status == NOT_STARTED and status_show:
            # keep track of start time/frame for later
            sfm_polygon_ordinate.frameNStart = frameN  # exact frame index
            sfm_polygon_ordinate.tStart = t  # local t and not account for scr refresh
            sfm_polygon_ordinate.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_polygon_ordinate, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_polygon_ordinate.status = STARTED
            sfm_polygon_ordinate.setAutoDraw(True)
        
        # if sfm_polygon_ordinate is active this frame...
        if sfm_polygon_ordinate.status == STARTED:
            # update params
            pass
        
        # *sfm_polygon_abscissa* updates
        
        # if sfm_polygon_abscissa is starting this frame...
        if sfm_polygon_abscissa.status == NOT_STARTED and pigeon_show:
            # keep track of start time/frame for later
            sfm_polygon_abscissa.frameNStart = frameN  # exact frame index
            sfm_polygon_abscissa.tStart = t  # local t and not account for scr refresh
            sfm_polygon_abscissa.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_polygon_abscissa, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_polygon_abscissa.status = STARTED
            sfm_polygon_abscissa.setAutoDraw(True)
        
        # if sfm_polygon_abscissa is active this frame...
        if sfm_polygon_abscissa.status == STARTED:
            # update params
            pass
        
        # *sfm_polygon_box* updates
        
        # if sfm_polygon_box is starting this frame...
        if sfm_polygon_box.status == NOT_STARTED and box_show:
            # keep track of start time/frame for later
            sfm_polygon_box.frameNStart = frameN  # exact frame index
            sfm_polygon_box.tStart = t  # local t and not account for scr refresh
            sfm_polygon_box.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_polygon_box, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_polygon_box.status = STARTED
            sfm_polygon_box.setAutoDraw(True)
        
        # if sfm_polygon_box is active this frame...
        if sfm_polygon_box.status == STARTED:
            # update params
            pass
        
        # *sfm_text_instruction* updates
        
        # if sfm_text_instruction is starting this frame...
        if sfm_text_instruction.status == NOT_STARTED and instruction_show:
            # keep track of start time/frame for later
            sfm_text_instruction.frameNStart = frameN  # exact frame index
            sfm_text_instruction.tStart = t  # local t and not account for scr refresh
            sfm_text_instruction.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sfm_text_instruction, 'tStartRefresh')  # time at next scr refresh
            # update status
            sfm_text_instruction.status = STARTED
            sfm_text_instruction.setAutoDraw(True)
        
        # if sfm_text_instruction is active this frame...
        if sfm_text_instruction.status == STARTED:
            # update params
            pass
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, win=win)
            return
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[routineTimer], 
                playbackComponents=[]
            )
            # skip the frame we paused on
            continue
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            show_final_message.forceEnded = routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in show_final_message.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "show_final_message" ---
    for thisComponent in show_final_message.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for show_final_message
    show_final_message.tStop = globalClock.getTime(format='float')
    show_final_message.tStopRefresh = tThisFlipGlobal
    thisExp.addData('show_final_message.stopped', show_final_message.tStop)
    # Run 'End Routine' code from sfm_code_update
    # Ick, but whatever. Makes sure that we don't ever
    # drop unexpectedly out of a move loop without 
    # resetting
    trial_start = True
    thisExp.nextEntry()
    # the Routine "show_final_message" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # mark experiment as finished
    endExperiment(thisExp, win=win)


def saveData(thisExp):
    """
    Save data from this experiment
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    filename = thisExp.dataFileName
    # these shouldn't be strictly necessary (should auto-save)
    thisExp.saveAsWideText(filename + '.csv', delim='auto')
    thisExp.saveAsPickle(filename)


def endExperiment(thisExp, win=None):
    """
    End this experiment, performing final shut down operations.
    
    This function does NOT close the window or end the Python process - use `quit` for this.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window for this experiment.
    """
    if win is not None:
        # remove autodraw from all current components
        win.clearAutoDraw()
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed
        win.flip()
    # return console logger level to WARNING
    logging.console.setLevel(logging.WARNING)
    # mark experiment handler as finished
    thisExp.status = FINISHED
    logging.flush()


def quit(thisExp, win=None, thisSession=None):
    """
    Fully quit, closing the window and ending the Python process.
    
    Parameters
    ==========
    win : psychopy.visual.Window
        Window to close.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    thisExp.abort()  # or data files will save again on exit
    # make sure everything is closed down
    if win is not None:
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed before quitting
        win.flip()
        win.close()
    logging.flush()
    if thisSession is not None:
        thisSession.stop()
    # terminate Python process
    core.quit()


# if running this experiment as a script...
if __name__ == '__main__':
    # call all functions in order
    expInfo = showExpInfoDlg(expInfo=expInfo)
    thisExp = setupData(expInfo=expInfo)
    logFile = setupLogging(filename=thisExp.dataFileName)
    win = setupWindow(expInfo=expInfo)
    setupDevices(expInfo=expInfo, thisExp=thisExp, win=win)
    run(
        expInfo=expInfo, 
        thisExp=thisExp, 
        win=win,
        globalClock='float'
    )
    saveData(thisExp=thisExp)
    quit(thisExp=thisExp, win=win)
