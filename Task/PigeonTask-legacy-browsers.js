/******************* 
 * Pigeontask Test *
 *******************/


// store info about the experiment session:
let expName = 'PigeonTask';  // from the Builder filename that created this script
let expInfo = {'participant': '', 'session': '001'};

// Start code blocks for 'Before Experiment'
// init psychoJS:
const psychoJS = new PsychoJS({
  debug: true
});

// open window:
psychoJS.openWindow({
  fullscr: false,
  color: new util.Color([0,0,0]),
  units: 'height',
  waitBlanking: true
});
// schedule the experiment:
psychoJS.schedule(psychoJS.gui.DlgFromDict({
  dictionary: expInfo,
  title: expName
}));

const flowScheduler = new Scheduler(psychoJS);
const dialogCancelScheduler = new Scheduler(psychoJS);
psychoJS.scheduleCondition(function() { return (psychoJS.gui.dialogComponent.button === 'OK'); }, flowScheduler, dialogCancelScheduler);

// flowScheduler gets run if the participants presses OK
flowScheduler.add(updateInfo); // add timeStamp
flowScheduler.add(experimentInit);
flowScheduler.add(consentRoutineBegin());
flowScheduler.add(consentRoutineEachFrame());
flowScheduler.add(consentRoutineEnd());
flowScheduler.add(experiment_setupRoutineBegin());
flowScheduler.add(experiment_setupRoutineEachFrame());
flowScheduler.add(experiment_setupRoutineEnd());
const general_instructions_loopLoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(general_instructions_loopLoopBegin(general_instructions_loopLoopScheduler));
flowScheduler.add(general_instructions_loopLoopScheduler);
flowScheduler.add(general_instructions_loopLoopEnd);
const simple_demo_block_loopLoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(simple_demo_block_loopLoopBegin(simple_demo_block_loopLoopScheduler));
flowScheduler.add(simple_demo_block_loopLoopScheduler);
flowScheduler.add(simple_demo_block_loopLoopEnd);
const status_instruction_loopLoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(status_instruction_loopLoopBegin(status_instruction_loopLoopScheduler));
flowScheduler.add(status_instruction_loopLoopScheduler);
flowScheduler.add(status_instruction_loopLoopEnd);
const experiment_loopLoopScheduler = new Scheduler(psychoJS);
flowScheduler.add(experiment_loopLoopBegin(experiment_loopLoopScheduler));
flowScheduler.add(experiment_loopLoopScheduler);
flowScheduler.add(experiment_loopLoopEnd);
flowScheduler.add(final_messageRoutineBegin());
flowScheduler.add(final_messageRoutineEachFrame());
flowScheduler.add(final_messageRoutineEnd());
flowScheduler.add(show_final_messageRoutineBegin());
flowScheduler.add(show_final_messageRoutineEachFrame());
flowScheduler.add(show_final_messageRoutineEnd());
flowScheduler.add(quitPsychoJS, '', true);

// quit if user presses Cancel in dialog box:
dialogCancelScheduler.add(quitPsychoJS, '', false);

psychoJS.start({
  expName: expName,
  expInfo: expInfo,
  resources: [
    {'name': 'Resources/PigeonTaskConditions.xlsx', 'path': 'Resources/PigeonTaskConditions.xlsx'},
    {'name': 'Resources/PigeonSNRSet01Conditions.xlsx', 'path': 'Resources/PigeonSNRSet01Conditions.xlsx'},
    {'name': 'Resources/PigeonSimpleDemoConditions.xlsx', 'path': 'Resources/PigeonSimpleDemoConditions.xlsx'},
    {'name': 'Resources/pigeon.png', 'path': 'Resources/pigeon.png'},
    {'name': 'Resources/Pigeon_predefined_commit_InstructionsConditions.xlsx', 'path': 'Resources/Pigeon_predefined_commit_InstructionsConditions.xlsx'},
    {'name': 'Resources/PigeonGeneralInstructionsConditions.xlsx', 'path': 'Resources/PigeonGeneralInstructionsConditions.xlsx'},
    {'name': 'Resources/PigeonStatusInstructionsConditions.xlsx', 'path': 'Resources/PigeonStatusInstructionsConditions.xlsx'},
    {'name': 'Resources/PigeonSNRSet02Conditions.xlsx', 'path': 'Resources/PigeonSNRSet02Conditions.xlsx'},
    {'name': 'Resources/seeds.png', 'path': 'Resources/seeds.png'}
  ]
});

psychoJS.experimentLogger.setLevel(core.Logger.ServerLevel.EXP);


var frameDur;
async function updateInfo() {
  expInfo['date'] = util.MonotonicClock.getDateStr();  // add a simple timestamp
  expInfo['expName'] = expName;
  expInfo['psychopyVersion'] = '2022.1.4';
  expInfo['OS'] = window.navigator.platform;

  psychoJS.experiment.dataFileName = (("." + "/") + `data/${expInfo["participant"]}_${expName}_${expInfo["date"]}`);

  // store frame rate of monitor if we can measure it successfully
  expInfo['frameRate'] = psychoJS.window.getActualFrameRate();
  if (typeof expInfo['frameRate'] !== 'undefined')
    frameDur = 1.0 / Math.round(expInfo['frameRate']);
  else
    frameDur = 1.0 / 60.0; // couldn't get a reliable measure so guess

  // add info from the URL:
  util.addInfoFromUrl(expInfo);
  psychoJS.setRedirectUrls('https://app.prolific.co/submissions/complete?cc=CKXJLI4I', '');

  return Scheduler.Event.NEXT;
}


var consentClock;
var consent_txt;
var consent_box;
var consent_mouse;
var consent_acknowledge;
var experiment_setupClock;
var still_going;
var choice;
var white;
var grey;
var yellow;
var green;
var lightgreen;
var black;
var red;
var ORIGIN_X;
var ORIGIN_Y;
var SEEDS_WIDTH;
var SEEDS_HEIGHT;
var EDGE_DISTANCE;
var INSTRUCTION_X;
var INSTRUCTION_Y;
var INSTRUCTION_COLOR;
var CHOICE_STRING;
var CONTINUE_STRING;
var PREDEFINED_CHOICE_STRING;
var PREDEFINED_CHOICE_COMMIT_STRING;
var PREDEFINED_NO_CHOICE_STRING;
var PREDEFINED_ABORT_INSTRUCTION;
var PREDEFINED_ABORT_STRING;
var CORRECT_STRING;
var ERROR_STRING;
var MESSAGE_X;
var MESSAGE_Y;
var instruction_show;
var instruction_string;
var message_show;
var message_string;
var steps_string;
var coins_string;
var LINE_HEIGHT;
var STATUS_BAR_HEIGHT;
var STATUS_BAR_WIDTH;
var STATUS_BAR_BASE_Y;
var STATUS_BAR_APEX_Y;
var STATUS_BAR_BRACKET_WIDTH;
var ABSCISSA_Y;
var ABSCISSA_WIDTH;
var ORDINATE_Y;
var ORDINATE_WIDTH;
var ORDINATE_HEIGHT;
var COINS_BAR_X;
var COINS_TEXT_X;
var COINS_TEXT_Y;
var COINS_LABEL_TEXT_X;
var COINS_LABEL_TEXT_Y;
var COINS_PROPERTIES_TEXT_X;
var COINS_PROPERTIES_TEXT_Y;
var STEPS_BAR_X;
var STEPS_TEXT_X;
var STEPS_TEXT_Y;
var STEPS_LABEL_TEXT_X;
var STEPS_LABEL_TEXT_Y;
var STEPS_PROPERTIES_TEXT_X;
var STEPS_PROPERTIES_TEXT_Y;
var steps_bar_height;
var steps_bar_center;
var coins_bar_height;
var coins_bar_center;
var PREDEFINED_BOUND_COLOR;
var PREDEFINED_BOUND_WIDTH;
var PREDEFINED_BOUND_HEIGHT;
var PREDEFINED_BOUND_Y;
var PREDEFINED_BOUND_SMALL_DELTA;
var PREDEFINED_BOUND_LARGE_DELTA;
var PREDEFINED_BOUND_X_DEFAULT;
var PREDEFINED_BOUND_MARKER_SZ;
var PREDEFINED_BOUND_MARKER_Y;
var PREDEFINED_BOUND_MARKER_COLOR;
var predefined_bound_x;
var predefined_bound_x_previous;
var committed_to_bound;
var predefined_bound_min_bonus;
var predefined_bound_max_bonus;
var min_steps_to_commit;
var seeds_show;
var pigeon_show;
var status_show;
var box_show;
var UPDATE_INTERVAL;
var PREDEFINED_UPDATE_INTERVAL;
var block_number;
var trial_number;
var trial_start;
var bonus_count;
var task_type;
var STEPS_MAX_PER_TRIAL;
var coins_gained_per_correct;
var coins_lost_per_error;
var coins_paid_to_start_trial;
var coins_count;
var coins_max;
var steps_lost_per_error;
var steps_taken_to_start_trial;
var steps_count;
var steps_max;
var step_mean_val;
var step_std_val;
var last_coins_count;
var last_steps_count;
var pigeon_flip;
var pigeon_true_x;
var pigeon_shown_x;
var pigeon_steps;
var move_pigeon;
var update_status;
var set_trial_end_string;
var set_block_end_string;
var general_instructions_setupClock;
var show_messageClock;
var sm_kbd_continue;
var box_width;
var box_height;
var box_x;
var box_y;
var sm_polygon_box;
var sm_polygon_abscissa;
var sm_polygon_ordinate;
var sm_polygon_petey_midpoint;
var sm_image_right_seeds;
var sm_image_left_seeds;
var sm_image_petey;
var sm_polygon_coin_bar;
var sm_polygon_steps_bar;
var sm_polygon_status_bar_apex;
var sm_polygon_status_bar_base;
var sm_text_coins;
var sm_text_coins_label;
var sm_text_coins_properties;
var sm_text_steps;
var sm_text_steps_label;
var sm_text_steps_properties;
var sm_text_message;
var sm_text_instruction;
var simple_demo_setupClock;
var simple_demoClock;
var sd_polygon_abscissa;
var sd_image_left_seeds;
var sd_image_right_seeds;
var sd_image_pigeon;
var status_instruction_setupClock;
var start_taskClock;
var setup_instructions_messageClock;
var online_task_runClock;
var otr_kbd_response;
var otr_polygon_abscissa;
var otr_polygon_ordinate;
var otr_polygon_petey_midpoint;
var otr_image_pigeon;
var otr_image_right_seeds;
var otr_image_left_seeds;
var otr_polygon_coins_bar;
var otr_polygon_steps_bar;
var otr_polygon_status_bar_apex;
var otr_polygon_status_bar_base;
var otr_text_coins;
var otr_text_coins_label;
var otr_text_coins_properties;
var otr_text_steps;
var otr_text_steps_label;
var otr_text_steps_properties;
var otr_text_instruction;
var predefined_task_setupClock;
var pts_kbd_update;
var pts_polygon_abscissa;
var pts_polygon_ordinate;
var pts_polygon_petey_midpoint;
var pts_polygon_left_bound;
var pts_polygon_right_bound;
var pts_polygon_stay_left;
var pts_polygon_stay_right;
var pts_image_pigeon;
var pts_image_right_seeds;
var pts_image_left_seeds;
var pts_polygon_coins_bar;
var pts_polygon_steps_bar;
var pts_polygon_status_bar_apex;
var pts_polygon_status_bar_base;
var pts_text_coins;
var pts_text_coins_label;
var pts_text_coins_properties;
var pts_text_steps;
var pts_text_steps_label;
var pts_text_steps_properties;
var pts_text_instruction;
var pts_text_hint;
var pts_text_hint_2;
var predefined_task_runClock;
var ptr_kbd_response;
var ptr_polygon_abscissa;
var ptr_polygon_ordinate;
var ptr_polygon_petey_midpoint;
var ptr_polygon_left_bound;
var ptr_polygon_right_bound;
var ptr_image_pigeon;
var ptr_image_right_seeds;
var ptr_image_left_seeds;
var ptr_polygon_coins_bar;
var ptr_polygon_steps_bar;
var ptr_polygon_status_bar_apex;
var ptr_polygon_status_bar_base;
var ptr_text_coins;
var ptr_text_coins_label;
var ptr_text_coins_properties;
var ptr_text_steps;
var ptr_text_steps_label;
var ptr_text_steps_properties;
var ptr_text_instruction;
var final_messageClock;
var show_final_messageClock;
var sfm_kbd_continue;
var sfm_text_message;
var sfm_text_steps_properties;
var sfm_text_steps_label;
var sfm_text_steps;
var sfm_text_coins_properties;
var sfm_text_coins_label;
var sfm_text_coins;
var sfm_polygon_status_bar_base;
var sfm_polygon_status_bar_apex;
var sfm_polygon_steps_bar;
var sfm_polygon_coin_bar;
var sfm_image_petey;
var sfm_image_left_seeds;
var sfm_image_right_seeds;
var sfm_polygon_petey_midpoint;
var sfm_polygon_ordinate;
var sfm_polygon_abscissa;
var sfm_polygon_box;
var sfm_text_instruction;
var globalClock;
var routineTimer;
async function experimentInit() {
  // Initialize components for Routine "consent"
  consentClock = new util.Clock();
  consent_txt = new visual.TextBox({
    win: psychoJS.window,
    name: 'consent_txt',
    text: 'This is an academic research project conducted through the University of Pennsylvania. In this game, a pigeon will wander from the center of the screen to one of two seed piles on the left or right. Your goal is either to guess the direction the pigeon will ultimately choose or set a cut-off for the pigeon’s movement to make a choice based on their current location.\n\nThe estimated total time is about 30 minutes.\n\nYou must be at least 18 years old to participate. Your participation in this research is voluntary, which means you can choose whether or not to participate without negative consequences. Your anonymity is assured: the researchers who have requested your participation will not receive any personal information about you except your previously provided Prolific demographic data such as gender, ethnicity, and age. The de-identified data may be stored and distributed for future research studies without additional informed consent.\n\nIf you have questions about this research, please contact Alice Dallstream by emailing adalls@pennmedicine.upenn.edu or Josh Gold by emailing jigold@pennmedicine.upenn.edu.\n\nIf you have questions, concerns, or complaints regarding your participation in this research study, or if you have any questions about your rights as a research subject and you cannot reach a member of the study team, you may contact the Office of Regulatory Affairs at the University of Pennsylvania by calling (215) 898-2614 or emailing irb@pobox.upenn.edu.\n',
    font: 'Open Sans',
    pos: [0, 0], letterHeight: 0.023,
    size: [0.75, 0.68],  units: undefined, 
    color: 'white', colorSpace: 'rgb',
    fillColor: [(- 0.1765), (- 0.1765), (- 0.1765)], borderColor: undefined,
    bold: false, italic: false,
    opacity: undefined,
    padding: 0.0,
    editable: false,
    multiline: true,
    anchor: 'center',
    depth: 0.0 
  });
  
  consent_box = new visual.Rect ({
    win: psychoJS.window, name: 'consent_box', 
    width: [0.05, 0.05][0], height: [0.05, 0.05][1],
    ori: 0.0, pos: [(- 0.4), (- 0.4)],
    lineWidth: 1.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color('white'),
    opacity: undefined, depth: -1, interpolate: true,
  });
  
  consent_mouse = new core.Mouse({
    win: psychoJS.window,
  });
  consent_mouse.mouseClock = new util.Clock();
  consent_acknowledge = new visual.TextBox({
    win: psychoJS.window,
    name: 'consent_acknowledge',
    text: 'I have read the information above and agree to take part in this study. I understand that I may withdraw my consent at any time before I complete the tasks.',
    font: 'Open Sans',
    pos: [(- 0.36), (- 0.42)], letterHeight: 0.023,
    size: [0.71, 0.08],  units: undefined, 
    color: 'white', colorSpace: 'rgb',
    fillColor: [(- 0.1765), (- 0.1765), (- 0.1765)], borderColor: undefined,
    bold: false, italic: false,
    opacity: undefined,
    padding: 0.0,
    editable: false,
    multiline: true,
    anchor: 'center-left',
    depth: -3.0 
  });
  
  // Initialize components for Routine "experiment_setup"
  experiment_setupClock = new util.Clock();
  still_going = true;
  choice = "";
  white = new util.Color([1, 1, 1]);
  grey = new util.Color([0.2, 0.2, 0.2]);
  yellow = new util.Color([1, 1, 0]);
  green = new util.Color([-1, 0, -1]);
  lightgreen = new util.Color([0, 1, 0]);
  black = new util.Color([-1, -1, -1]);
  red = new util.Color([1, 0, 0]);
  ORIGIN_X = 0;
  ORIGIN_Y = 0;
  SEEDS_WIDTH = 0.1;
  SEEDS_HEIGHT = 0.1;
  EDGE_DISTANCE = 0.75;
  INSTRUCTION_X = ORIGIN_X;
  INSTRUCTION_Y = 0.45;
  INSTRUCTION_COLOR = yellow;
  CHOICE_STRING = "Press Left or Right arrow to choose";
  CONTINUE_STRING = "Press Spacebar or any arrow to continue";
  PREDEFINED_CHOICE_STRING = "Press Left/Up for small/large increases, Right/Down for small/large decreases,\nSpacebar to start the trial.";
  PREDEFINED_CHOICE_COMMIT_STRING = "Press Left/Up for small/large increases, Right/Down for small/large decreases,\nc to commit or Spacebar to start the trial.";
  PREDEFINED_NO_CHOICE_STRING = "Petey failed to reach the bound\nConsider moving it closer.";
  PREDEFINED_ABORT_INSTRUCTION = "Press Spacebar or any arrow to abort trial";
  PREDEFINED_ABORT_STRING = "Aborted early";
  CORRECT_STRING = "Correct";
  ERROR_STRING = "Error";
  MESSAGE_X = ORIGIN_X;
  MESSAGE_Y = 0.3;
  instruction_show = true;
  instruction_string = CONTINUE_STRING;
  message_show = true;
  message_string = "";
  steps_string = " ";
  coins_string = " ";
  LINE_HEIGHT = 0.1;
  STATUS_BAR_HEIGHT = 0.2;
  STATUS_BAR_WIDTH = 0.1;
  STATUS_BAR_BASE_Y = (ORIGIN_Y - 0.4);
  STATUS_BAR_APEX_Y = (STATUS_BAR_BASE_Y + STATUS_BAR_HEIGHT);
  STATUS_BAR_BRACKET_WIDTH = 0.35;
  ABSCISSA_Y = (ORIGIN_Y - 0.05);
  ABSCISSA_WIDTH = (EDGE_DISTANCE * 2);
  ORDINATE_Y = (ORIGIN_Y + 0.05);
  ORDINATE_WIDTH = 0.2;
  ORDINATE_HEIGHT = 0.01;
  COINS_BAR_X = (- 0.1);
  COINS_TEXT_X = (- 0.1);
  COINS_TEXT_Y = (- 0.37);
  COINS_LABEL_TEXT_X = (- 0.28);
  COINS_LABEL_TEXT_Y = (- 0.2);
  COINS_PROPERTIES_TEXT_X = (- 0.55);
  COINS_PROPERTIES_TEXT_Y = (- 0.3);
  STEPS_BAR_X = 0.1;
  STEPS_TEXT_X = 0.1;
  STEPS_TEXT_Y = (- 0.37);
  STEPS_LABEL_TEXT_X = 0.3;
  STEPS_LABEL_TEXT_Y = (- 0.2);
  STEPS_PROPERTIES_TEXT_X = 0.55;
  STEPS_PROPERTIES_TEXT_Y = (- 0.3);
  steps_bar_height = 0;
  steps_bar_center = 0;
  coins_bar_height = 0;
  coins_bar_center = 0;
  PREDEFINED_BOUND_COLOR = green;
  PREDEFINED_BOUND_WIDTH = 0.2;
  PREDEFINED_BOUND_HEIGHT = 0.01;
  PREDEFINED_BOUND_Y = (ORIGIN_Y + 0.05);
  PREDEFINED_BOUND_SMALL_DELTA = 0.003;
  PREDEFINED_BOUND_LARGE_DELTA = 0.05;
  PREDEFINED_BOUND_X_DEFAULT = 0.05;
  PREDEFINED_BOUND_MARKER_SZ = 0.01;
  PREDEFINED_BOUND_MARKER_Y = 0.16;
  PREDEFINED_BOUND_MARKER_COLOR = lightgreen;
  predefined_bound_x = PREDEFINED_BOUND_X_DEFAULT;
  predefined_bound_x_previous = predefined_bound_x;
  committed_to_bound = false;
  predefined_bound_min_bonus = 0;
  predefined_bound_max_bonus = 0;
  min_steps_to_commit = 100;
  seeds_show = false;
  pigeon_show = false;
  status_show = false;
  box_show = false;
  UPDATE_INTERVAL = 0.2;
  PREDEFINED_UPDATE_INTERVAL = 0.1;
  block_number = 0;
  trial_number = 0;
  trial_start = true;
  bonus_count = 0;
  task_type = "online";
  STEPS_MAX_PER_TRIAL = 50;
  coins_gained_per_correct = 1;
  coins_lost_per_error = 1;
  coins_paid_to_start_trial = 0;
  coins_count = 2;
  coins_max = 8;
  steps_lost_per_error = 0;
  steps_taken_to_start_trial = 1;
  steps_count = 0;
  steps_max = 200;
  step_mean_val = 0.02;
  step_std_val = 0.1;
  last_coins_count = 0;
  last_steps_count = 0;
  pigeon_flip = 0.1;
  pigeon_true_x = 0;
  pigeon_shown_x = pigeon_true_x;
  pigeon_steps = [];
  function move_pigeon_function(step_mean, step_std) {
      var new_step, r, u, v;
      for (var i, _pj_c = 0, _pj_a = util.range(10), _pj_b = _pj_a.length; (_pj_c < _pj_b); _pj_c += 1) {
          i = _pj_a[_pj_c];
          u = ((2 * Math.random()) - 1);
          v = ((2 * Math.random()) - 1);
          r = ((u * u) + (v * v));
          if (((r !== 0) && (r < 1))) {
              break;
          }
      }
      new_step = (step_mean + ((step_std * u) * Math.sqrt((((- 2) * Math.log(r)) / r))));
      if ((new_step > 0)) {
          pigeon_flip = (- 0.1);
      } else {
          pigeon_flip = 0.1;
      }
      pigeon_true_x = (pigeon_true_x + new_step);
      if ((pigeon_true_x >= EDGE_DISTANCE)) {
          pigeon_shown_x = EDGE_DISTANCE;
      } else {
          if ((pigeon_true_x <= (- EDGE_DISTANCE))) {
              pigeon_shown_x = (- EDGE_DISTANCE);
          } else {
              pigeon_shown_x = pigeon_true_x;
          }
      }
  }
  function update_status_function(update_text, update_coins_bar, update_steps_bar) {
      if (update_text) {
          coins_string = ((((((("Coins gained per correct choice = " + coins_gained_per_correct.toString()) + "\nCoins lost per error choice = ") + coins_lost_per_error.toString()) + "\nCoins paid to start each trial = ") + coins_paid_to_start_trial.toString()) + "\nTarget number of coins for bonus = ") + coins_max.toString());
          steps_string = ((((("Total steps = " + steps_max.toString()) + "\nSteps lost per error choice = ") + steps_lost_per_error.toString()) + "\nSteps taken to start each trial = ") + steps_taken_to_start_trial.toString());
      }
      if (update_coins_bar) {
          if ((coins_count <= coins_max)) {
              coins_bar_height = ((coins_count / coins_max) * STATUS_BAR_HEIGHT);
          } else {
              coins_bar_height = (STATUS_BAR_HEIGHT + 0.02);
          }
          coins_bar_center = (STATUS_BAR_BASE_Y + (coins_bar_height / 2));
      }
      if (update_steps_bar) {
          steps_bar_height = (((steps_max - steps_count) / steps_max) * STATUS_BAR_HEIGHT);
          steps_bar_center = (STATUS_BAR_BASE_Y + (steps_bar_height / 2));
      }
  }
  function set_trial_end_string_function(feedback_string) {
      var coins_change, coins_feedback_string, steps_change;
      coins_change = (coins_count - last_coins_count);
      steps_change = (steps_count - last_steps_count);
      if ((coins_change > 0)) {
          coins_feedback_string = (coins_change.toString() + " coins gained");
      } else {
          if ((coins_change < 0)) {
              coins_feedback_string = (coins_change.toString() + " coins lost");
          } else {
              coins_feedback_string = "no coins gained or lost";
          }
      }
      message_string = (((((feedback_string + "\n\n") + coins_feedback_string) + "\n\n") + steps_change.toString()) + " steps used this trial");
      last_coins_count = coins_count;
      last_steps_count = steps_count;
  }
  function set_block_end_string_function() {
      if (((((task_type === "predefined_commit") && (predefined_bound_x >= predefined_bound_min_bonus)) && (predefined_bound_x <= predefined_bound_max_bonus)) || ((task_type !== "predefined_commit") && (coins_count >= coins_max)))) {
          message_string = "Block finished.\n\nYou earned a bonus. Great job!";
          bonus_count += 1;
          coins_bar_height = STATUS_BAR_HEIGHT;
          coins_bar_center = (STATUS_BAR_BASE_Y + (coins_bar_height / 2));
          coins_count = "Bonus!";
      } else {
          message_string = "Block finished.\n\nGreat job!";
          coins_bar_height = 0;
          coins_bar_center = STATUS_BAR_BASE_Y;
          coins_count = "No bonus";
      }
  }
  move_pigeon = move_pigeon_function;
  update_status = update_status_function;
  set_trial_end_string = set_trial_end_string_function;
  set_block_end_string = set_block_end_string_function;
  update_status(true, true, true);
  
  // Initialize components for Routine "general_instructions_setup"
  general_instructions_setupClock = new util.Clock();
  // Initialize components for Routine "show_message"
  show_messageClock = new util.Clock();
  sm_kbd_continue = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  box_width = 0;
  box_height = 0;
  box_x = 0;
  box_y = 0;
  
  sm_polygon_box = new visual.Rect ({
    win: psychoJS.window, name: 'sm_polygon_box', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 3.0, lineColor: new util.Color([1.0, (- 1.0), (- 1.0)]),
    fillColor: new util.Color([0.0, 0.0, 0.0]),
    opacity: undefined, depth: -2, interpolate: true,
  });
  
  sm_polygon_abscissa = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sm_polygon_abscissa', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -3, interpolate: true,
  });
  
  sm_polygon_ordinate = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sm_polygon_ordinate', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -4, interpolate: true,
  });
  
  sm_polygon_petey_midpoint = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sm_polygon_petey_midpoint', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([1.0, (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -5, interpolate: true,
  });
  
  sm_image_right_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'sm_image_right_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [EDGE_DISTANCE, ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -6.0 
  });
  sm_image_left_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'sm_image_left_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [(- EDGE_DISTANCE), ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -7.0 
  });
  sm_image_petey = new visual.ImageStim({
    win : psychoJS.window,
    name : 'sm_image_petey', units : undefined, 
    image : 'Resources/pigeon.png', mask : undefined,
    ori : 0.0, pos : [0, 0], size : 1.0,
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -8.0 
  });
  sm_polygon_coin_bar = new visual.Rect ({
    win: psychoJS.window, name: 'sm_polygon_coin_bar', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 1.0, lineColor: new util.Color([1.0, 1.0, 1.0]),
    fillColor: new util.Color([1.0, 1.0, 1.0]),
    opacity: undefined, depth: -9, interpolate: true,
  });
  
  sm_polygon_steps_bar = new visual.Rect ({
    win: psychoJS.window, name: 'sm_polygon_steps_bar', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 1.0, lineColor: new util.Color([1.0, 1.0, 1.0]),
    fillColor: new util.Color([1.0, 1.0, 1.0]),
    opacity: undefined, depth: -10, interpolate: true,
  });
  
  sm_polygon_status_bar_apex = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sm_polygon_status_bar_apex', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -11, interpolate: true,
  });
  
  sm_polygon_status_bar_base = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sm_polygon_status_bar_base', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -12, interpolate: true,
  });
  
  sm_text_coins = new visual.TextStim({
    win: psychoJS.window,
    name: 'sm_text_coins',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_TEXT_X, COINS_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),  opacity: undefined,
    depth: -13.0 
  });
  
  sm_text_coins_label = new visual.TextStim({
    win: psychoJS.window,
    name: 'sm_text_coins_label',
    text: 'Coins earned',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_LABEL_TEXT_X, COINS_LABEL_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),  opacity: undefined,
    depth: -14.0 
  });
  
  sm_text_coins_properties = new visual.TextStim({
    win: psychoJS.window,
    name: 'sm_text_coins_properties',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_PROPERTIES_TEXT_X, COINS_PROPERTIES_TEXT_Y], height: 0.025,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([1.0, 1.0, 1.0]),  opacity: undefined,
    depth: -15.0 
  });
  
  sm_text_steps = new visual.TextStim({
    win: psychoJS.window,
    name: 'sm_text_steps',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_TEXT_X, STEPS_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),  opacity: undefined,
    depth: -16.0 
  });
  
  sm_text_steps_label = new visual.TextStim({
    win: psychoJS.window,
    name: 'sm_text_steps_label',
    text: 'Steps remaining',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),  opacity: undefined,
    depth: -17.0 
  });
  
  sm_text_steps_properties = new visual.TextStim({
    win: psychoJS.window,
    name: 'sm_text_steps_properties',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_PROPERTIES_TEXT_X, STEPS_PROPERTIES_TEXT_Y], height: 0.025,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([1.0, 1.0, 1.0]),  opacity: undefined,
    depth: -18.0 
  });
  
  sm_text_message = new visual.TextStim({
    win: psychoJS.window,
    name: 'sm_text_message',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [MESSAGE_X, MESSAGE_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([1.0, 1.0, 1.0]),  opacity: undefined,
    depth: -19.0 
  });
  
  sm_text_instruction = new visual.TextStim({
    win: psychoJS.window,
    name: 'sm_text_instruction',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [INSTRUCTION_X, INSTRUCTION_Y], height: 0.02,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color(INSTRUCTION_COLOR),  opacity: undefined,
    depth: -20.0 
  });
  
  // Initialize components for Routine "simple_demo_setup"
  simple_demo_setupClock = new util.Clock();
  // Initialize components for Routine "simple_demo"
  simple_demoClock = new util.Clock();
  sd_polygon_abscissa = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sd_polygon_abscissa', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -1, interpolate: true,
  });
  
  sd_image_left_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'sd_image_left_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [(- EDGE_DISTANCE), ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -2.0 
  });
  sd_image_right_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'sd_image_right_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [EDGE_DISTANCE, ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -3.0 
  });
  sd_image_pigeon = new visual.ImageStim({
    win : psychoJS.window,
    name : 'sd_image_pigeon', units : undefined, 
    image : 'Resources/pigeon.png', mask : undefined,
    ori : 0.0, pos : [0, 0], size : 1.0,
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -4.0 
  });
  // Initialize components for Routine "status_instruction_setup"
  status_instruction_setupClock = new util.Clock();
  // Initialize components for Routine "start_task"
  start_taskClock = new util.Clock();
  // Initialize components for Routine "setup_instructions_message"
  setup_instructions_messageClock = new util.Clock();
  // Initialize components for Routine "online_task_run"
  online_task_runClock = new util.Clock();
  otr_kbd_response = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  otr_polygon_abscissa = new visual.ShapeStim ({
    win: psychoJS.window, name: 'otr_polygon_abscissa', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -2, interpolate: true,
  });
  
  otr_polygon_ordinate = new visual.ShapeStim ({
    win: psychoJS.window, name: 'otr_polygon_ordinate', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -3, interpolate: true,
  });
  
  otr_polygon_petey_midpoint = new visual.ShapeStim ({
    win: psychoJS.window, name: 'otr_polygon_petey_midpoint', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color('red'),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -4, interpolate: true,
  });
  
  otr_image_pigeon = new visual.ImageStim({
    win : psychoJS.window,
    name : 'otr_image_pigeon', units : undefined, 
    image : 'Resources/pigeon.png', mask : undefined,
    ori : 0.0, pos : [0, 0], size : 1.0,
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -5.0 
  });
  otr_image_right_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'otr_image_right_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [EDGE_DISTANCE, ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -6.0 
  });
  otr_image_left_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'otr_image_left_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [(- EDGE_DISTANCE), ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -7.0 
  });
  otr_polygon_coins_bar = new visual.Rect ({
    win: psychoJS.window, name: 'otr_polygon_coins_bar', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 1.0, lineColor: new util.Color('white'),
    fillColor: new util.Color('white'),
    opacity: undefined, depth: -8, interpolate: true,
  });
  
  otr_polygon_steps_bar = new visual.Rect ({
    win: psychoJS.window, name: 'otr_polygon_steps_bar', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 1.0, lineColor: new util.Color('white'),
    fillColor: new util.Color('white'),
    opacity: undefined, depth: -9, interpolate: true,
  });
  
  otr_polygon_status_bar_apex = new visual.ShapeStim ({
    win: psychoJS.window, name: 'otr_polygon_status_bar_apex', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color('black'),
    fillColor: new util.Color('black'),
    opacity: undefined, depth: -10, interpolate: true,
  });
  
  otr_polygon_status_bar_base = new visual.ShapeStim ({
    win: psychoJS.window, name: 'otr_polygon_status_bar_base', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -11, interpolate: true,
  });
  
  otr_text_coins = new visual.TextStim({
    win: psychoJS.window,
    name: 'otr_text_coins',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_TEXT_X, COINS_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -12.0 
  });
  
  otr_text_coins_label = new visual.TextStim({
    win: psychoJS.window,
    name: 'otr_text_coins_label',
    text: 'Coins earned',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_LABEL_TEXT_X, COINS_LABEL_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -13.0 
  });
  
  otr_text_coins_properties = new visual.TextStim({
    win: psychoJS.window,
    name: 'otr_text_coins_properties',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_PROPERTIES_TEXT_X, COINS_PROPERTIES_TEXT_Y], height: 0.025,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('white'),  opacity: undefined,
    depth: -14.0 
  });
  
  otr_text_steps = new visual.TextStim({
    win: psychoJS.window,
    name: 'otr_text_steps',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_TEXT_X, STEPS_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -15.0 
  });
  
  otr_text_steps_label = new visual.TextStim({
    win: psychoJS.window,
    name: 'otr_text_steps_label',
    text: 'Steps remaining',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -16.0 
  });
  
  otr_text_steps_properties = new visual.TextStim({
    win: psychoJS.window,
    name: 'otr_text_steps_properties',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_PROPERTIES_TEXT_X, STEPS_PROPERTIES_TEXT_Y], height: 0.025,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('white'),  opacity: undefined,
    depth: -17.0 
  });
  
  otr_text_instruction = new visual.TextStim({
    win: psychoJS.window,
    name: 'otr_text_instruction',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [INSTRUCTION_X, INSTRUCTION_Y], height: 0.02,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color(INSTRUCTION_COLOR),  opacity: undefined,
    depth: -18.0 
  });
  
  // Initialize components for Routine "predefined_task_setup"
  predefined_task_setupClock = new util.Clock();
  pts_kbd_update = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  pts_polygon_abscissa = new visual.ShapeStim ({
    win: psychoJS.window, name: 'pts_polygon_abscissa', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -2, interpolate: true,
  });
  
  pts_polygon_ordinate = new visual.ShapeStim ({
    win: psychoJS.window, name: 'pts_polygon_ordinate', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -3, interpolate: true,
  });
  
  pts_polygon_petey_midpoint = new visual.ShapeStim ({
    win: psychoJS.window, name: 'pts_polygon_petey_midpoint', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color('red'),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -4, interpolate: true,
  });
  
  pts_polygon_left_bound = new visual.ShapeStim ({
    win: psychoJS.window, name: 'pts_polygon_left_bound', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 3.0, lineColor: new util.Color('white'),
    fillColor: new util.Color([0.0, 0.0, 0.0]),
    opacity: undefined, depth: -5, interpolate: true,
  });
  
  pts_polygon_right_bound = new visual.ShapeStim ({
    win: psychoJS.window, name: 'pts_polygon_right_bound', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 3.0, lineColor: new util.Color(PREDEFINED_BOUND_COLOR),
    fillColor: new util.Color([0.0, 0.0, 0.0]),
    opacity: undefined, depth: -6, interpolate: true,
  });
  
  pts_polygon_stay_left = new visual.ShapeStim ({
    win: psychoJS.window, name: 'pts_polygon_stay_left', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color(PREDEFINED_BOUND_MARKER_COLOR),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -7, interpolate: true,
  });
  
  pts_polygon_stay_right = new visual.ShapeStim ({
    win: psychoJS.window, name: 'pts_polygon_stay_right', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color(PREDEFINED_BOUND_MARKER_COLOR),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -8, interpolate: true,
  });
  
  pts_image_pigeon = new visual.ImageStim({
    win : psychoJS.window,
    name : 'pts_image_pigeon', units : undefined, 
    image : 'Resources/pigeon.png', mask : undefined,
    ori : 0.0, pos : [0, 0], size : 1.0,
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -9.0 
  });
  pts_image_right_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'pts_image_right_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [EDGE_DISTANCE, ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -10.0 
  });
  pts_image_left_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'pts_image_left_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [(- EDGE_DISTANCE), ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -11.0 
  });
  pts_polygon_coins_bar = new visual.Rect ({
    win: psychoJS.window, name: 'pts_polygon_coins_bar', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 1.0, lineColor: new util.Color('white'),
    fillColor: new util.Color('white'),
    opacity: undefined, depth: -12, interpolate: true,
  });
  
  pts_polygon_steps_bar = new visual.Rect ({
    win: psychoJS.window, name: 'pts_polygon_steps_bar', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 1.0, lineColor: new util.Color('white'),
    fillColor: new util.Color('white'),
    opacity: undefined, depth: -13, interpolate: true,
  });
  
  pts_polygon_status_bar_apex = new visual.ShapeStim ({
    win: psychoJS.window, name: 'pts_polygon_status_bar_apex', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color('black'),
    fillColor: new util.Color('black'),
    opacity: undefined, depth: -14, interpolate: true,
  });
  
  pts_polygon_status_bar_base = new visual.ShapeStim ({
    win: psychoJS.window, name: 'pts_polygon_status_bar_base', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -15, interpolate: true,
  });
  
  pts_text_coins = new visual.TextStim({
    win: psychoJS.window,
    name: 'pts_text_coins',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [0, 0], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -16.0 
  });
  
  pts_text_coins_label = new visual.TextStim({
    win: psychoJS.window,
    name: 'pts_text_coins_label',
    text: 'Coins earned',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_LABEL_TEXT_X, COINS_LABEL_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -17.0 
  });
  
  pts_text_coins_properties = new visual.TextStim({
    win: psychoJS.window,
    name: 'pts_text_coins_properties',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_PROPERTIES_TEXT_X, COINS_PROPERTIES_TEXT_Y], height: 0.025,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('white'),  opacity: undefined,
    depth: -18.0 
  });
  
  pts_text_steps = new visual.TextStim({
    win: psychoJS.window,
    name: 'pts_text_steps',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_TEXT_X, STEPS_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -19.0 
  });
  
  pts_text_steps_label = new visual.TextStim({
    win: psychoJS.window,
    name: 'pts_text_steps_label',
    text: 'Steps remaining',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -20.0 
  });
  
  pts_text_steps_properties = new visual.TextStim({
    win: psychoJS.window,
    name: 'pts_text_steps_properties',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_PROPERTIES_TEXT_X, STEPS_PROPERTIES_TEXT_Y], height: 0.025,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('white'),  opacity: undefined,
    depth: -21.0 
  });
  
  pts_text_instruction = new visual.TextStim({
    win: psychoJS.window,
    name: 'pts_text_instruction',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [INSTRUCTION_X, INSTRUCTION_Y], height: 0.02,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color(INSTRUCTION_COLOR),  opacity: undefined,
    depth: -22.0 
  });
  
  pts_text_hint = new visual.TextStim({
    win: psychoJS.window,
    name: 'pts_text_hint',
    text: '-> move closer to save steps <-\n<- move farther to increase accuracy ->',
    font: 'Open Sans',
    units: undefined, 
    pos: [0, 0], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('white'),  opacity: undefined,
    depth: -23.0 
  });
  
  pts_text_hint_2 = new visual.TextStim({
    win: psychoJS.window,
    name: 'pts_text_hint_2',
    text: 'The light green tics indicate your last choice.',
    font: 'Open Sans',
    units: undefined, 
    pos: [0, 0], height: 0.02,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([0.1294, 0.8667, 0.1294]),  opacity: undefined,
    depth: -24.0 
  });
  
  // Initialize components for Routine "predefined_task_run"
  predefined_task_runClock = new util.Clock();
  ptr_kbd_response = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  ptr_polygon_abscissa = new visual.ShapeStim ({
    win: psychoJS.window, name: 'ptr_polygon_abscissa', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -2, interpolate: true,
  });
  
  ptr_polygon_ordinate = new visual.ShapeStim ({
    win: psychoJS.window, name: 'ptr_polygon_ordinate', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -3, interpolate: true,
  });
  
  ptr_polygon_petey_midpoint = new visual.ShapeStim ({
    win: psychoJS.window, name: 'ptr_polygon_petey_midpoint', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color('red'),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -4, interpolate: true,
  });
  
  ptr_polygon_left_bound = new visual.ShapeStim ({
    win: psychoJS.window, name: 'ptr_polygon_left_bound', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 3.0, lineColor: new util.Color('white'),
    fillColor: new util.Color([0.0, 0.0, 0.0]),
    opacity: undefined, depth: -5, interpolate: true,
  });
  
  ptr_polygon_right_bound = new visual.ShapeStim ({
    win: psychoJS.window, name: 'ptr_polygon_right_bound', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 3.0, lineColor: new util.Color('white'),
    fillColor: new util.Color([0.0, 0.0, 0.0]),
    opacity: undefined, depth: -6, interpolate: true,
  });
  
  ptr_image_pigeon = new visual.ImageStim({
    win : psychoJS.window,
    name : 'ptr_image_pigeon', units : undefined, 
    image : 'Resources/pigeon.png', mask : undefined,
    ori : 0.0, pos : [0, 0], size : 1.0,
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -7.0 
  });
  ptr_image_right_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'ptr_image_right_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [EDGE_DISTANCE, ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -8.0 
  });
  ptr_image_left_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'ptr_image_left_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [(- EDGE_DISTANCE), ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -9.0 
  });
  ptr_polygon_coins_bar = new visual.Rect ({
    win: psychoJS.window, name: 'ptr_polygon_coins_bar', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 1.0, lineColor: new util.Color('white'),
    fillColor: new util.Color('white'),
    opacity: undefined, depth: -10, interpolate: true,
  });
  
  ptr_polygon_steps_bar = new visual.Rect ({
    win: psychoJS.window, name: 'ptr_polygon_steps_bar', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 1.0, lineColor: new util.Color('white'),
    fillColor: new util.Color('white'),
    opacity: undefined, depth: -11, interpolate: true,
  });
  
  ptr_polygon_status_bar_apex = new visual.ShapeStim ({
    win: psychoJS.window, name: 'ptr_polygon_status_bar_apex', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color('black'),
    fillColor: new util.Color('black'),
    opacity: undefined, depth: -12, interpolate: true,
  });
  
  ptr_polygon_status_bar_base = new visual.ShapeStim ({
    win: psychoJS.window, name: 'ptr_polygon_status_bar_base', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -13, interpolate: true,
  });
  
  ptr_text_coins = new visual.TextStim({
    win: psychoJS.window,
    name: 'ptr_text_coins',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_TEXT_X, COINS_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -14.0 
  });
  
  ptr_text_coins_label = new visual.TextStim({
    win: psychoJS.window,
    name: 'ptr_text_coins_label',
    text: 'Coins earned',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_LABEL_TEXT_X, COINS_LABEL_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -15.0 
  });
  
  ptr_text_coins_properties = new visual.TextStim({
    win: psychoJS.window,
    name: 'ptr_text_coins_properties',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_PROPERTIES_TEXT_X, COINS_PROPERTIES_TEXT_Y], height: 0.025,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('white'),  opacity: undefined,
    depth: -16.0 
  });
  
  ptr_text_steps = new visual.TextStim({
    win: psychoJS.window,
    name: 'ptr_text_steps',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_TEXT_X, STEPS_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -17.0 
  });
  
  ptr_text_steps_label = new visual.TextStim({
    win: psychoJS.window,
    name: 'ptr_text_steps_label',
    text: 'Steps remaining',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('black'),  opacity: undefined,
    depth: -18.0 
  });
  
  ptr_text_steps_properties = new visual.TextStim({
    win: psychoJS.window,
    name: 'ptr_text_steps_properties',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_PROPERTIES_TEXT_X, STEPS_PROPERTIES_TEXT_Y], height: 0.025,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color('white'),  opacity: undefined,
    depth: -19.0 
  });
  
  ptr_text_instruction = new visual.TextStim({
    win: psychoJS.window,
    name: 'ptr_text_instruction',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [INSTRUCTION_X, INSTRUCTION_Y], height: 0.02,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color(INSTRUCTION_COLOR),  opacity: undefined,
    depth: -20.0 
  });
  
  // Initialize components for Routine "final_message"
  final_messageClock = new util.Clock();
  // Initialize components for Routine "show_final_message"
  show_final_messageClock = new util.Clock();
  sfm_kbd_continue = new core.Keyboard({psychoJS: psychoJS, clock: new util.Clock(), waitForStart: true});
  
  box_width = 0;
  box_height = 0;
  box_x = 0;
  box_y = 0;
  
  sfm_text_message = new visual.TextStim({
    win: psychoJS.window,
    name: 'sfm_text_message',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [MESSAGE_X, MESSAGE_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([1.0, 1.0, 1.0]),  opacity: undefined,
    depth: -2.0 
  });
  
  sfm_text_steps_properties = new visual.TextStim({
    win: psychoJS.window,
    name: 'sfm_text_steps_properties',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_PROPERTIES_TEXT_X, STEPS_PROPERTIES_TEXT_Y], height: 0.025,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([1.0, 1.0, 1.0]),  opacity: undefined,
    depth: -3.0 
  });
  
  sfm_text_steps_label = new visual.TextStim({
    win: psychoJS.window,
    name: 'sfm_text_steps_label',
    text: 'Steps remaining',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_LABEL_TEXT_X, STEPS_LABEL_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),  opacity: undefined,
    depth: -4.0 
  });
  
  sfm_text_steps = new visual.TextStim({
    win: psychoJS.window,
    name: 'sfm_text_steps',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [STEPS_TEXT_X, STEPS_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),  opacity: undefined,
    depth: -5.0 
  });
  
  sfm_text_coins_properties = new visual.TextStim({
    win: psychoJS.window,
    name: 'sfm_text_coins_properties',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_PROPERTIES_TEXT_X, COINS_PROPERTIES_TEXT_Y], height: 0.025,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([1.0, 1.0, 1.0]),  opacity: undefined,
    depth: -6.0 
  });
  
  sfm_text_coins_label = new visual.TextStim({
    win: psychoJS.window,
    name: 'sfm_text_coins_label',
    text: 'Coins earned',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_LABEL_TEXT_X, COINS_LABEL_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),  opacity: undefined,
    depth: -7.0 
  });
  
  sfm_text_coins = new visual.TextStim({
    win: psychoJS.window,
    name: 'sfm_text_coins',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [COINS_TEXT_X, COINS_TEXT_Y], height: 0.03,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),  opacity: undefined,
    depth: -8.0 
  });
  
  sfm_polygon_status_bar_base = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sfm_polygon_status_bar_base', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -9, interpolate: true,
  });
  
  sfm_polygon_status_bar_apex = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sfm_polygon_status_bar_apex', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -10, interpolate: true,
  });
  
  sfm_polygon_steps_bar = new visual.Rect ({
    win: psychoJS.window, name: 'sfm_polygon_steps_bar', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 1.0, lineColor: new util.Color([1.0, 1.0, 1.0]),
    fillColor: new util.Color([1.0, 1.0, 1.0]),
    opacity: undefined, depth: -11, interpolate: true,
  });
  
  sfm_polygon_coin_bar = new visual.Rect ({
    win: psychoJS.window, name: 'sfm_polygon_coin_bar', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 1.0, lineColor: new util.Color([1.0, 1.0, 1.0]),
    fillColor: new util.Color([1.0, 1.0, 1.0]),
    opacity: undefined, depth: -12, interpolate: true,
  });
  
  sfm_image_petey = new visual.ImageStim({
    win : psychoJS.window,
    name : 'sfm_image_petey', units : undefined, 
    image : 'Resources/pigeon.png', mask : undefined,
    ori : 0.0, pos : [0, 0], size : 1.0,
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -13.0 
  });
  sfm_image_left_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'sfm_image_left_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [(- EDGE_DISTANCE), ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -14.0 
  });
  sfm_image_right_seeds = new visual.ImageStim({
    win : psychoJS.window,
    name : 'sfm_image_right_seeds', units : undefined, 
    image : 'Resources/seeds.png', mask : undefined,
    ori : 0.0, pos : [EDGE_DISTANCE, ORIGIN_Y], size : [SEEDS_WIDTH, SEEDS_HEIGHT],
    color : new util.Color([1,1,1]), opacity : undefined,
    flipHoriz : false, flipVert : false,
    texRes : 128.0, interpolate : true, depth : -15.0 
  });
  sfm_polygon_petey_midpoint = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sfm_polygon_petey_midpoint', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([1.0, (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -16, interpolate: true,
  });
  
  sfm_polygon_ordinate = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sfm_polygon_ordinate', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 90.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -17, interpolate: true,
  });
  
  sfm_polygon_abscissa = new visual.ShapeStim ({
    win: psychoJS.window, name: 'sfm_polygon_abscissa', 
    vertices: [[-[1.0, 1.0][0]/2.0, 0], [+[1.0, 1.0][0]/2.0, 0]],
    ori: 0.0, pos: [0, 0],
    lineWidth: 2.0, lineColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    fillColor: new util.Color([(- 1.0), (- 1.0), (- 1.0)]),
    opacity: undefined, depth: -18, interpolate: true,
  });
  
  sfm_polygon_box = new visual.Rect ({
    win: psychoJS.window, name: 'sfm_polygon_box', 
    width: [1.0, 1.0][0], height: [1.0, 1.0][1],
    ori: 0.0, pos: [0, 0],
    lineWidth: 3.0, lineColor: new util.Color([1.0, (- 1.0), (- 1.0)]),
    fillColor: new util.Color([0.0, 0.0, 0.0]),
    opacity: undefined, depth: -19, interpolate: true,
  });
  
  sfm_text_instruction = new visual.TextStim({
    win: psychoJS.window,
    name: 'sfm_text_instruction',
    text: '',
    font: 'Open Sans',
    units: undefined, 
    pos: [INSTRUCTION_X, INSTRUCTION_Y], height: 0.02,  wrapWidth: undefined, ori: 0.0,
    color: new util.Color(INSTRUCTION_COLOR),  opacity: undefined,
    depth: -20.0 
  });
  
  // Create some handy timers
  globalClock = new util.Clock();  // to track the time since experiment started
  routineTimer = new util.CountdownTimer();  // to track time remaining of each (non-slip) routine
  
  return Scheduler.Event.NEXT;
}


var t;
var frameN;
var continueRoutine;
var gotValidClick;
var consentComponents;
function consentRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'consent'-------
    t = 0;
    consentClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    // setup some python lists for storing info about the consent_mouse
    consent_mouse.clicked_name = [];
    gotValidClick = false; // until a click is received
    // keep track of which components have finished
    consentComponents = [];
    consentComponents.push(consent_txt);
    consentComponents.push(consent_box);
    consentComponents.push(consent_mouse);
    consentComponents.push(consent_acknowledge);
    
    consentComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


var prevButtonState;
var _mouseButtons;
function consentRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'consent'-------
    // get current time
    t = consentClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *consent_txt* updates
    if (t >= 0.0 && consent_txt.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      consent_txt.tStart = t;  // (not accounting for frame time here)
      consent_txt.frameNStart = frameN;  // exact frame index
      
      consent_txt.setAutoDraw(true);
    }

    
    // *consent_box* updates
    if (t >= 0.0 && consent_box.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      consent_box.tStart = t;  // (not accounting for frame time here)
      consent_box.frameNStart = frameN;  // exact frame index
      
      consent_box.setAutoDraw(true);
    }

    // *consent_mouse* updates
    if (t >= 0.0 && consent_mouse.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      consent_mouse.tStart = t;  // (not accounting for frame time here)
      consent_mouse.frameNStart = frameN;  // exact frame index
      
      consent_mouse.status = PsychoJS.Status.STARTED;
      consent_mouse.mouseClock.reset();
      prevButtonState = consent_mouse.getPressed();  // if button is down already this ISN'T a new click
      }
    if (consent_mouse.status === PsychoJS.Status.STARTED) {  // only update if started and not finished!
      _mouseButtons = consent_mouse.getPressed();
      if (!_mouseButtons.every( (e,i,) => (e == prevButtonState[i]) )) { // button state changed?
        prevButtonState = _mouseButtons;
        if (_mouseButtons.reduce( (e, acc) => (e+acc) ) > 0) { // state changed to a new click
          // check if the mouse was inside our 'clickable' objects
          gotValidClick = false;
          for (const obj of [consent_box]) {
            if (obj.contains(consent_mouse)) {
              gotValidClick = true;
              consent_mouse.clicked_name.push(obj.name)
            }
          }
          if (gotValidClick === true) { // abort routine on response
            continueRoutine = false;
          }
        }
      }
    }
    
    // *consent_acknowledge* updates
    if (t >= 0.0 && consent_acknowledge.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      consent_acknowledge.tStart = t;  // (not accounting for frame time here)
      consent_acknowledge.frameNStart = frameN;  // exact frame index
      
      consent_acknowledge.setAutoDraw(true);
    }

    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    consentComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


var _mouseXYs;
function consentRoutineEnd() {
  return async function () {
    //------Ending Routine 'consent'-------
    consentComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    // store data for psychoJS.experiment (ExperimentHandler)
    _mouseXYs = consent_mouse.getPos();
    _mouseButtons = consent_mouse.getPressed();
    psychoJS.experiment.addData('consent_mouse.x', _mouseXYs[0]);
    psychoJS.experiment.addData('consent_mouse.y', _mouseXYs[1]);
    psychoJS.experiment.addData('consent_mouse.leftButton', _mouseButtons[0]);
    psychoJS.experiment.addData('consent_mouse.midButton', _mouseButtons[1]);
    psychoJS.experiment.addData('consent_mouse.rightButton', _mouseButtons[2]);
    if (consent_mouse.clicked_name.length > 0) {
      psychoJS.experiment.addData('consent_mouse.clicked_name', consent_mouse.clicked_name[0]);}
    // the Routine "consent" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var experiment_setupComponents;
function experiment_setupRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'experiment_setup'-------
    t = 0;
    experiment_setupClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    // keep track of which components have finished
    experiment_setupComponents = [];
    
    experiment_setupComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function experiment_setupRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'experiment_setup'-------
    // get current time
    t = experiment_setupClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    experiment_setupComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function experiment_setupRoutineEnd() {
  return async function () {
    //------Ending Routine 'experiment_setup'-------
    experiment_setupComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    // the Routine "experiment_setup" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var general_instructions_loop;
var currentLoop;
function general_instructions_loopLoopBegin(general_instructions_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    general_instructions_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 1, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: 'Resources/PigeonGeneralInstructionsConditions.xlsx',
      seed: undefined, name: 'general_instructions_loop'
    });
    psychoJS.experiment.addLoop(general_instructions_loop); // add the loop to the experiment
    currentLoop = general_instructions_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    general_instructions_loop.forEach(function() {
      const snapshot = general_instructions_loop.getSnapshot();
    
      general_instructions_loopLoopScheduler.add(importConditions(snapshot));
      general_instructions_loopLoopScheduler.add(general_instructions_setupRoutineBegin(snapshot));
      general_instructions_loopLoopScheduler.add(general_instructions_setupRoutineEachFrame());
      general_instructions_loopLoopScheduler.add(general_instructions_setupRoutineEnd());
      general_instructions_loopLoopScheduler.add(show_messageRoutineBegin(snapshot));
      general_instructions_loopLoopScheduler.add(show_messageRoutineEachFrame());
      general_instructions_loopLoopScheduler.add(show_messageRoutineEnd());
      general_instructions_loopLoopScheduler.add(endLoopIteration(general_instructions_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


async function general_instructions_loopLoopEnd() {
  psychoJS.experiment.removeLoop(general_instructions_loop);

  return Scheduler.Event.NEXT;
}


var simple_demo_block_loop;
function simple_demo_block_loopLoopBegin(simple_demo_block_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    simple_demo_block_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 1, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: 'Resources/PigeonSimpleDemoConditions.xlsx',
      seed: undefined, name: 'simple_demo_block_loop'
    });
    psychoJS.experiment.addLoop(simple_demo_block_loop); // add the loop to the experiment
    currentLoop = simple_demo_block_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    simple_demo_block_loop.forEach(function() {
      const snapshot = simple_demo_block_loop.getSnapshot();
    
      simple_demo_block_loopLoopScheduler.add(importConditions(snapshot));
      simple_demo_block_loopLoopScheduler.add(simple_demo_setupRoutineBegin(snapshot));
      simple_demo_block_loopLoopScheduler.add(simple_demo_setupRoutineEachFrame());
      simple_demo_block_loopLoopScheduler.add(simple_demo_setupRoutineEnd());
      simple_demo_block_loopLoopScheduler.add(show_messageRoutineBegin(snapshot));
      simple_demo_block_loopLoopScheduler.add(show_messageRoutineEachFrame());
      simple_demo_block_loopLoopScheduler.add(show_messageRoutineEnd());
      const simple_demo_trial_loopLoopScheduler = new Scheduler(psychoJS);
      simple_demo_block_loopLoopScheduler.add(simple_demo_trial_loopLoopBegin(simple_demo_trial_loopLoopScheduler, snapshot));
      simple_demo_block_loopLoopScheduler.add(simple_demo_trial_loopLoopScheduler);
      simple_demo_block_loopLoopScheduler.add(simple_demo_trial_loopLoopEnd);
      simple_demo_block_loopLoopScheduler.add(endLoopIteration(simple_demo_block_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


var simple_demo_trial_loop;
function simple_demo_trial_loopLoopBegin(simple_demo_trial_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    simple_demo_trial_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 1000, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: undefined,
      seed: undefined, name: 'simple_demo_trial_loop'
    });
    psychoJS.experiment.addLoop(simple_demo_trial_loop); // add the loop to the experiment
    currentLoop = simple_demo_trial_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    simple_demo_trial_loop.forEach(function() {
      const snapshot = simple_demo_trial_loop.getSnapshot();
    
      simple_demo_trial_loopLoopScheduler.add(importConditions(snapshot));
      simple_demo_trial_loopLoopScheduler.add(simple_demoRoutineBegin(snapshot));
      simple_demo_trial_loopLoopScheduler.add(simple_demoRoutineEachFrame());
      simple_demo_trial_loopLoopScheduler.add(simple_demoRoutineEnd());
      simple_demo_trial_loopLoopScheduler.add(endLoopIteration(simple_demo_trial_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


async function simple_demo_trial_loopLoopEnd() {
  psychoJS.experiment.removeLoop(simple_demo_trial_loop);

  return Scheduler.Event.NEXT;
}


async function simple_demo_block_loopLoopEnd() {
  psychoJS.experiment.removeLoop(simple_demo_block_loop);

  return Scheduler.Event.NEXT;
}


var status_instruction_loop;
function status_instruction_loopLoopBegin(status_instruction_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    status_instruction_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 1, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: 'Resources/PigeonStatusInstructionsConditions.xlsx',
      seed: undefined, name: 'status_instruction_loop'
    });
    psychoJS.experiment.addLoop(status_instruction_loop); // add the loop to the experiment
    currentLoop = status_instruction_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    status_instruction_loop.forEach(function() {
      const snapshot = status_instruction_loop.getSnapshot();
    
      status_instruction_loopLoopScheduler.add(importConditions(snapshot));
      status_instruction_loopLoopScheduler.add(status_instruction_setupRoutineBegin(snapshot));
      status_instruction_loopLoopScheduler.add(status_instruction_setupRoutineEachFrame());
      status_instruction_loopLoopScheduler.add(status_instruction_setupRoutineEnd());
      status_instruction_loopLoopScheduler.add(show_messageRoutineBegin(snapshot));
      status_instruction_loopLoopScheduler.add(show_messageRoutineEachFrame());
      status_instruction_loopLoopScheduler.add(show_messageRoutineEnd());
      status_instruction_loopLoopScheduler.add(endLoopIteration(status_instruction_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


async function status_instruction_loopLoopEnd() {
  psychoJS.experiment.removeLoop(status_instruction_loop);

  return Scheduler.Event.NEXT;
}


var experiment_loop;
function experiment_loopLoopBegin(experiment_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    experiment_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 1, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: 'Resources/PigeonTaskConditions.xlsx',
      seed: undefined, name: 'experiment_loop'
    });
    psychoJS.experiment.addLoop(experiment_loop); // add the loop to the experiment
    currentLoop = experiment_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    experiment_loop.forEach(function() {
      const snapshot = experiment_loop.getSnapshot();
    
      experiment_loopLoopScheduler.add(importConditions(snapshot));
      experiment_loopLoopScheduler.add(start_taskRoutineBegin(snapshot));
      experiment_loopLoopScheduler.add(start_taskRoutineEachFrame());
      experiment_loopLoopScheduler.add(start_taskRoutineEnd());
      experiment_loopLoopScheduler.add(show_messageRoutineBegin(snapshot));
      experiment_loopLoopScheduler.add(show_messageRoutineEachFrame());
      experiment_loopLoopScheduler.add(show_messageRoutineEnd());
      const instructions_loopLoopScheduler = new Scheduler(psychoJS);
      experiment_loopLoopScheduler.add(instructions_loopLoopBegin(instructions_loopLoopScheduler, snapshot));
      experiment_loopLoopScheduler.add(instructions_loopLoopScheduler);
      experiment_loopLoopScheduler.add(instructions_loopLoopEnd);
      const online_trial_loopLoopScheduler = new Scheduler(psychoJS);
      experiment_loopLoopScheduler.add(online_trial_loopLoopBegin(online_trial_loopLoopScheduler, snapshot));
      experiment_loopLoopScheduler.add(online_trial_loopLoopScheduler);
      experiment_loopLoopScheduler.add(online_trial_loopLoopEnd);
      const predefined_trial_loopLoopScheduler = new Scheduler(psychoJS);
      experiment_loopLoopScheduler.add(predefined_trial_loopLoopBegin(predefined_trial_loopLoopScheduler, snapshot));
      experiment_loopLoopScheduler.add(predefined_trial_loopLoopScheduler);
      experiment_loopLoopScheduler.add(predefined_trial_loopLoopEnd);
      experiment_loopLoopScheduler.add(endLoopIteration(experiment_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


var instructions_loop;
function instructions_loopLoopBegin(instructions_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    instructions_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 1, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: instructions_conditions_file,
      seed: undefined, name: 'instructions_loop'
    });
    psychoJS.experiment.addLoop(instructions_loop); // add the loop to the experiment
    currentLoop = instructions_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    instructions_loop.forEach(function() {
      const snapshot = instructions_loop.getSnapshot();
    
      instructions_loopLoopScheduler.add(importConditions(snapshot));
      instructions_loopLoopScheduler.add(setup_instructions_messageRoutineBegin(snapshot));
      instructions_loopLoopScheduler.add(setup_instructions_messageRoutineEachFrame());
      instructions_loopLoopScheduler.add(setup_instructions_messageRoutineEnd());
      instructions_loopLoopScheduler.add(show_messageRoutineBegin(snapshot));
      instructions_loopLoopScheduler.add(show_messageRoutineEachFrame());
      instructions_loopLoopScheduler.add(show_messageRoutineEnd());
      instructions_loopLoopScheduler.add(endLoopIteration(instructions_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


async function instructions_loopLoopEnd() {
  psychoJS.experiment.removeLoop(instructions_loop);

  return Scheduler.Event.NEXT;
}


var online_trial_loop;
function online_trial_loopLoopBegin(online_trial_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    online_trial_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: online_trials_count, method: TrialHandler.Method.FULLRANDOM,
      extraInfo: expInfo, originPath: undefined,
      trialList: snr_conditions_file,
      seed: undefined, name: 'online_trial_loop'
    });
    psychoJS.experiment.addLoop(online_trial_loop); // add the loop to the experiment
    currentLoop = online_trial_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    online_trial_loop.forEach(function() {
      const snapshot = online_trial_loop.getSnapshot();
    
      online_trial_loopLoopScheduler.add(importConditions(snapshot));
      const online_move_loopLoopScheduler = new Scheduler(psychoJS);
      online_trial_loopLoopScheduler.add(online_move_loopLoopBegin(online_move_loopLoopScheduler, snapshot));
      online_trial_loopLoopScheduler.add(online_move_loopLoopScheduler);
      online_trial_loopLoopScheduler.add(online_move_loopLoopEnd);
      online_trial_loopLoopScheduler.add(show_messageRoutineBegin(snapshot));
      online_trial_loopLoopScheduler.add(show_messageRoutineEachFrame());
      online_trial_loopLoopScheduler.add(show_messageRoutineEnd());
      online_trial_loopLoopScheduler.add(endLoopIteration(online_trial_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


var online_move_loop;
function online_move_loopLoopBegin(online_move_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    online_move_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: STEPS_MAX_PER_TRIAL, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: undefined,
      seed: undefined, name: 'online_move_loop'
    });
    psychoJS.experiment.addLoop(online_move_loop); // add the loop to the experiment
    currentLoop = online_move_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    online_move_loop.forEach(function() {
      const snapshot = online_move_loop.getSnapshot();
    
      online_move_loopLoopScheduler.add(importConditions(snapshot));
      online_move_loopLoopScheduler.add(online_task_runRoutineBegin(snapshot));
      online_move_loopLoopScheduler.add(online_task_runRoutineEachFrame());
      online_move_loopLoopScheduler.add(online_task_runRoutineEnd());
      online_move_loopLoopScheduler.add(endLoopIteration(online_move_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


async function online_move_loopLoopEnd() {
  psychoJS.experiment.removeLoop(online_move_loop);

  return Scheduler.Event.NEXT;
}


async function online_trial_loopLoopEnd() {
  psychoJS.experiment.removeLoop(online_trial_loop);

  return Scheduler.Event.NEXT;
}


var predefined_trial_loop;
function predefined_trial_loopLoopBegin(predefined_trial_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    predefined_trial_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: predefined_trials_count, method: TrialHandler.Method.RANDOM,
      extraInfo: expInfo, originPath: undefined,
      trialList: snr_conditions_file,
      seed: undefined, name: 'predefined_trial_loop'
    });
    psychoJS.experiment.addLoop(predefined_trial_loop); // add the loop to the experiment
    currentLoop = predefined_trial_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    predefined_trial_loop.forEach(function() {
      const snapshot = predefined_trial_loop.getSnapshot();
    
      predefined_trial_loopLoopScheduler.add(importConditions(snapshot));
      const predefined_setup_loopLoopScheduler = new Scheduler(psychoJS);
      predefined_trial_loopLoopScheduler.add(predefined_setup_loopLoopBegin(predefined_setup_loopLoopScheduler, snapshot));
      predefined_trial_loopLoopScheduler.add(predefined_setup_loopLoopScheduler);
      predefined_trial_loopLoopScheduler.add(predefined_setup_loopLoopEnd);
      const predefined_move_loopLoopScheduler = new Scheduler(psychoJS);
      predefined_trial_loopLoopScheduler.add(predefined_move_loopLoopBegin(predefined_move_loopLoopScheduler, snapshot));
      predefined_trial_loopLoopScheduler.add(predefined_move_loopLoopScheduler);
      predefined_trial_loopLoopScheduler.add(predefined_move_loopLoopEnd);
      predefined_trial_loopLoopScheduler.add(show_messageRoutineBegin(snapshot));
      predefined_trial_loopLoopScheduler.add(show_messageRoutineEachFrame());
      predefined_trial_loopLoopScheduler.add(show_messageRoutineEnd());
      predefined_trial_loopLoopScheduler.add(endLoopIteration(predefined_trial_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


var predefined_setup_loop;
function predefined_setup_loopLoopBegin(predefined_setup_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    predefined_setup_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: 1000, method: TrialHandler.Method.SEQUENTIAL,
      extraInfo: expInfo, originPath: undefined,
      trialList: undefined,
      seed: undefined, name: 'predefined_setup_loop'
    });
    psychoJS.experiment.addLoop(predefined_setup_loop); // add the loop to the experiment
    currentLoop = predefined_setup_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    predefined_setup_loop.forEach(function() {
      const snapshot = predefined_setup_loop.getSnapshot();
    
      predefined_setup_loopLoopScheduler.add(importConditions(snapshot));
      predefined_setup_loopLoopScheduler.add(predefined_task_setupRoutineBegin(snapshot));
      predefined_setup_loopLoopScheduler.add(predefined_task_setupRoutineEachFrame());
      predefined_setup_loopLoopScheduler.add(predefined_task_setupRoutineEnd());
      predefined_setup_loopLoopScheduler.add(endLoopIteration(predefined_setup_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


async function predefined_setup_loopLoopEnd() {
  psychoJS.experiment.removeLoop(predefined_setup_loop);

  return Scheduler.Event.NEXT;
}


var predefined_move_loop;
function predefined_move_loopLoopBegin(predefined_move_loopLoopScheduler, snapshot) {
  return async function() {
    TrialHandler.fromSnapshot(snapshot); // update internal variables (.thisN etc) of the loop
    
    // set up handler to look after randomisation of conditions etc
    predefined_move_loop = new TrialHandler({
      psychoJS: psychoJS,
      nReps: STEPS_MAX_PER_TRIAL, method: TrialHandler.Method.RANDOM,
      extraInfo: expInfo, originPath: undefined,
      trialList: undefined,
      seed: undefined, name: 'predefined_move_loop'
    });
    psychoJS.experiment.addLoop(predefined_move_loop); // add the loop to the experiment
    currentLoop = predefined_move_loop;  // we're now the current loop
    
    // Schedule all the trials in the trialList:
    predefined_move_loop.forEach(function() {
      const snapshot = predefined_move_loop.getSnapshot();
    
      predefined_move_loopLoopScheduler.add(importConditions(snapshot));
      predefined_move_loopLoopScheduler.add(predefined_task_runRoutineBegin(snapshot));
      predefined_move_loopLoopScheduler.add(predefined_task_runRoutineEachFrame());
      predefined_move_loopLoopScheduler.add(predefined_task_runRoutineEnd());
      predefined_move_loopLoopScheduler.add(endLoopIteration(predefined_move_loopLoopScheduler, snapshot));
    });
    
    return Scheduler.Event.NEXT;
  }
}


async function predefined_move_loopLoopEnd() {
  psychoJS.experiment.removeLoop(predefined_move_loop);

  return Scheduler.Event.NEXT;
}


async function predefined_trial_loopLoopEnd() {
  psychoJS.experiment.removeLoop(predefined_trial_loop);

  return Scheduler.Event.NEXT;
}


async function experiment_loopLoopEnd() {
  psychoJS.experiment.removeLoop(experiment_loop);

  return Scheduler.Event.NEXT;
}


var general_instructions_setupComponents;
function general_instructions_setupRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'general_instructions_setup'-------
    t = 0;
    general_instructions_setupClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    status_show = false;
    box_show = false;
    pigeon_show = gic_pigeon_show;
    seeds_show = gic_seeds_show;
    message_string = gic_message;
    
    // keep track of which components have finished
    general_instructions_setupComponents = [];
    
    general_instructions_setupComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function general_instructions_setupRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'general_instructions_setup'-------
    // get current time
    t = general_instructions_setupClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    general_instructions_setupComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function general_instructions_setupRoutineEnd() {
  return async function () {
    //------Ending Routine 'general_instructions_setup'-------
    general_instructions_setupComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    // the Routine "general_instructions_setup" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _sm_kbd_continue_allKeys;
var show_messageComponents;
function show_messageRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'show_message'-------
    t = 0;
    show_messageClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    sm_kbd_continue.keys = undefined;
    sm_kbd_continue.rt = undefined;
    _sm_kbd_continue_allKeys = [];
    sm_kbd_continue.clearEvents();
    
    sm_polygon_box.setPos([box_x, box_y]);
    sm_polygon_box.setSize([box_width, box_height]);
    sm_polygon_abscissa.setPos([ORIGIN_X, ABSCISSA_Y]);
    sm_polygon_abscissa.setSize([ABSCISSA_WIDTH, LINE_HEIGHT]);
    sm_polygon_ordinate.setPos([ORIGIN_X, ORDINATE_Y]);
    sm_polygon_ordinate.setSize([ORDINATE_WIDTH, ORDINATE_HEIGHT]);
    sm_polygon_petey_midpoint.setPos([pigeon_shown_x, ORDINATE_Y]);
    sm_polygon_petey_midpoint.setSize([(ORDINATE_WIDTH * 0.9), ORDINATE_HEIGHT]);
    sm_image_petey.setPos([pigeon_shown_x, ORIGIN_Y]);
    sm_image_petey.setSize([pigeon_flip, 0.1]);
    sm_polygon_coin_bar.setPos([COINS_BAR_X, coins_bar_center]);
    sm_polygon_coin_bar.setSize([STATUS_BAR_WIDTH, coins_bar_height]);
    sm_polygon_steps_bar.setPos([STEPS_BAR_X, steps_bar_center]);
    sm_polygon_steps_bar.setSize([STATUS_BAR_WIDTH, steps_bar_height]);
    sm_polygon_status_bar_apex.setPos([ORIGIN_X, STATUS_BAR_APEX_Y]);
    sm_polygon_status_bar_apex.setSize([STATUS_BAR_BRACKET_WIDTH, LINE_HEIGHT]);
    sm_polygon_status_bar_base.setPos([ORIGIN_X, STATUS_BAR_BASE_Y]);
    sm_polygon_status_bar_base.setSize([STATUS_BAR_BRACKET_WIDTH, LINE_HEIGHT]);
    sm_text_coins.setText(coins_count);
    sm_text_coins_properties.setText(coins_string);
    sm_text_steps.setText((steps_max - steps_count));
    sm_text_steps_properties.setText(steps_string);
    sm_text_message.setText(message_string);
    sm_text_instruction.setText(CONTINUE_STRING);
    // keep track of which components have finished
    show_messageComponents = [];
    show_messageComponents.push(sm_kbd_continue);
    show_messageComponents.push(sm_polygon_box);
    show_messageComponents.push(sm_polygon_abscissa);
    show_messageComponents.push(sm_polygon_ordinate);
    show_messageComponents.push(sm_polygon_petey_midpoint);
    show_messageComponents.push(sm_image_right_seeds);
    show_messageComponents.push(sm_image_left_seeds);
    show_messageComponents.push(sm_image_petey);
    show_messageComponents.push(sm_polygon_coin_bar);
    show_messageComponents.push(sm_polygon_steps_bar);
    show_messageComponents.push(sm_polygon_status_bar_apex);
    show_messageComponents.push(sm_polygon_status_bar_base);
    show_messageComponents.push(sm_text_coins);
    show_messageComponents.push(sm_text_coins_label);
    show_messageComponents.push(sm_text_coins_properties);
    show_messageComponents.push(sm_text_steps);
    show_messageComponents.push(sm_text_steps_label);
    show_messageComponents.push(sm_text_steps_properties);
    show_messageComponents.push(sm_text_message);
    show_messageComponents.push(sm_text_instruction);
    
    show_messageComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function show_messageRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'show_message'-------
    // get current time
    t = show_messageClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *sm_kbd_continue* updates
    if (t >= 0.0 && sm_kbd_continue.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_kbd_continue.tStart = t;  // (not accounting for frame time here)
      sm_kbd_continue.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { sm_kbd_continue.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { sm_kbd_continue.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { sm_kbd_continue.clearEvents(); });
    }

    if (sm_kbd_continue.status === PsychoJS.Status.STARTED) {
      let theseKeys = sm_kbd_continue.getKeys({keyList: ['space', 'up', 'down', 'left', 'right'], waitRelease: false});
      _sm_kbd_continue_allKeys = _sm_kbd_continue_allKeys.concat(theseKeys);
      if (_sm_kbd_continue_allKeys.length > 0) {
        sm_kbd_continue.keys = _sm_kbd_continue_allKeys[_sm_kbd_continue_allKeys.length - 1].name;  // just the last key pressed
        sm_kbd_continue.rt = _sm_kbd_continue_allKeys[_sm_kbd_continue_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *sm_polygon_box* updates
    if ((box_show) && sm_polygon_box.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_polygon_box.tStart = t;  // (not accounting for frame time here)
      sm_polygon_box.frameNStart = frameN;  // exact frame index
      
      sm_polygon_box.setAutoDraw(true);
    }

    
    // *sm_polygon_abscissa* updates
    if ((pigeon_show) && sm_polygon_abscissa.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_polygon_abscissa.tStart = t;  // (not accounting for frame time here)
      sm_polygon_abscissa.frameNStart = frameN;  // exact frame index
      
      sm_polygon_abscissa.setAutoDraw(true);
    }

    
    // *sm_polygon_ordinate* updates
    if ((status_show) && sm_polygon_ordinate.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_polygon_ordinate.tStart = t;  // (not accounting for frame time here)
      sm_polygon_ordinate.frameNStart = frameN;  // exact frame index
      
      sm_polygon_ordinate.setAutoDraw(true);
    }

    
    // *sm_polygon_petey_midpoint* updates
    if ((status_show) && sm_polygon_petey_midpoint.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_polygon_petey_midpoint.tStart = t;  // (not accounting for frame time here)
      sm_polygon_petey_midpoint.frameNStart = frameN;  // exact frame index
      
      sm_polygon_petey_midpoint.setAutoDraw(true);
    }

    
    // *sm_image_right_seeds* updates
    if ((seeds_show) && sm_image_right_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_image_right_seeds.tStart = t;  // (not accounting for frame time here)
      sm_image_right_seeds.frameNStart = frameN;  // exact frame index
      
      sm_image_right_seeds.setAutoDraw(true);
    }

    
    // *sm_image_left_seeds* updates
    if ((seeds_show) && sm_image_left_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_image_left_seeds.tStart = t;  // (not accounting for frame time here)
      sm_image_left_seeds.frameNStart = frameN;  // exact frame index
      
      sm_image_left_seeds.setAutoDraw(true);
    }

    
    // *sm_image_petey* updates
    if ((pigeon_show) && sm_image_petey.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_image_petey.tStart = t;  // (not accounting for frame time here)
      sm_image_petey.frameNStart = frameN;  // exact frame index
      
      sm_image_petey.setAutoDraw(true);
    }

    
    // *sm_polygon_coin_bar* updates
    if ((status_show) && sm_polygon_coin_bar.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_polygon_coin_bar.tStart = t;  // (not accounting for frame time here)
      sm_polygon_coin_bar.frameNStart = frameN;  // exact frame index
      
      sm_polygon_coin_bar.setAutoDraw(true);
    }

    
    // *sm_polygon_steps_bar* updates
    if ((status_show) && sm_polygon_steps_bar.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_polygon_steps_bar.tStart = t;  // (not accounting for frame time here)
      sm_polygon_steps_bar.frameNStart = frameN;  // exact frame index
      
      sm_polygon_steps_bar.setAutoDraw(true);
    }

    
    // *sm_polygon_status_bar_apex* updates
    if ((status_show) && sm_polygon_status_bar_apex.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_polygon_status_bar_apex.tStart = t;  // (not accounting for frame time here)
      sm_polygon_status_bar_apex.frameNStart = frameN;  // exact frame index
      
      sm_polygon_status_bar_apex.setAutoDraw(true);
    }

    
    // *sm_polygon_status_bar_base* updates
    if ((status_show) && sm_polygon_status_bar_base.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_polygon_status_bar_base.tStart = t;  // (not accounting for frame time here)
      sm_polygon_status_bar_base.frameNStart = frameN;  // exact frame index
      
      sm_polygon_status_bar_base.setAutoDraw(true);
    }

    
    // *sm_text_coins* updates
    if ((status_show) && sm_text_coins.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_text_coins.tStart = t;  // (not accounting for frame time here)
      sm_text_coins.frameNStart = frameN;  // exact frame index
      
      sm_text_coins.setAutoDraw(true);
    }

    
    // *sm_text_coins_label* updates
    if ((status_show) && sm_text_coins_label.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_text_coins_label.tStart = t;  // (not accounting for frame time here)
      sm_text_coins_label.frameNStart = frameN;  // exact frame index
      
      sm_text_coins_label.setAutoDraw(true);
    }

    
    // *sm_text_coins_properties* updates
    if ((status_show) && sm_text_coins_properties.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_text_coins_properties.tStart = t;  // (not accounting for frame time here)
      sm_text_coins_properties.frameNStart = frameN;  // exact frame index
      
      sm_text_coins_properties.setAutoDraw(true);
    }

    
    // *sm_text_steps* updates
    if ((status_show) && sm_text_steps.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_text_steps.tStart = t;  // (not accounting for frame time here)
      sm_text_steps.frameNStart = frameN;  // exact frame index
      
      sm_text_steps.setAutoDraw(true);
    }

    
    // *sm_text_steps_label* updates
    if ((status_show) && sm_text_steps_label.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_text_steps_label.tStart = t;  // (not accounting for frame time here)
      sm_text_steps_label.frameNStart = frameN;  // exact frame index
      
      sm_text_steps_label.setAutoDraw(true);
    }

    
    // *sm_text_steps_properties* updates
    if ((status_show) && sm_text_steps_properties.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_text_steps_properties.tStart = t;  // (not accounting for frame time here)
      sm_text_steps_properties.frameNStart = frameN;  // exact frame index
      
      sm_text_steps_properties.setAutoDraw(true);
    }

    
    // *sm_text_message* updates
    if ((message_show) && sm_text_message.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_text_message.tStart = t;  // (not accounting for frame time here)
      sm_text_message.frameNStart = frameN;  // exact frame index
      
      sm_text_message.setAutoDraw(true);
    }

    
    // *sm_text_instruction* updates
    if ((instruction_show) && sm_text_instruction.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sm_text_instruction.tStart = t;  // (not accounting for frame time here)
      sm_text_instruction.frameNStart = frameN;  // exact frame index
      
      sm_text_instruction.setAutoDraw(true);
    }

    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    show_messageComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function show_messageRoutineEnd() {
  return async function () {
    //------Ending Routine 'show_message'-------
    show_messageComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    sm_kbd_continue.stop();
    trial_start = true;
    
    // the Routine "show_message" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var demo_start;
var step_mean;
var step_std;
var simple_demo_setupComponents;
function simple_demo_setupRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'simple_demo_setup'-------
    t = 0;
    simple_demo_setupClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    seeds_show = true;
    pigeon_show = true;
    status_show = false;
    box_show = false;
    demo_start = true;
    message_string = sdc_message;
    step_mean = sdc_step_mean;
    step_std = sdc_step_std;
    
    // keep track of which components have finished
    simple_demo_setupComponents = [];
    
    simple_demo_setupComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function simple_demo_setupRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'simple_demo_setup'-------
    // get current time
    t = simple_demo_setupClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    simple_demo_setupComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function simple_demo_setupRoutineEnd() {
  return async function () {
    //------Ending Routine 'simple_demo_setup'-------
    simple_demo_setupComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    // the Routine "simple_demo_setup" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var simple_demoComponents;
function simple_demoRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'simple_demo'-------
    t = 0;
    simple_demoClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    if (demo_start) {
        pigeon_flip = 0.1;
        pigeon_true_x = 0;
        pigeon_shown_x = pigeon_true_x;
        if ((util.randint(0, 2) === 0)) {
            step_mean_val = ((- 1) * Math.abs(step_mean));
        } else {
            step_mean_val = Math.abs(step_mean);
        }
        step_std_val = step_mean;
        demo_start = false;
    } else {
        if (((pigeon_true_x > (- EDGE_DISTANCE)) && (pigeon_true_x < EDGE_DISTANCE))) {
            move_pigeon(step_mean_val, step_std_val);
        } else {
            simple_demo_trial_loop.finished = true;
        }
    }
    
    sd_polygon_abscissa.setPos([ORIGIN_X, ABSCISSA_Y]);
    sd_polygon_abscissa.setSize([(EDGE_DISTANCE * 2), LINE_HEIGHT]);
    sd_image_pigeon.setPos([pigeon_shown_x, ORIGIN_Y]);
    sd_image_pigeon.setSize([pigeon_flip, 0.1]);
    // keep track of which components have finished
    simple_demoComponents = [];
    simple_demoComponents.push(sd_polygon_abscissa);
    simple_demoComponents.push(sd_image_left_seeds);
    simple_demoComponents.push(sd_image_right_seeds);
    simple_demoComponents.push(sd_image_pigeon);
    
    simple_demoComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


var frameRemains;
function simple_demoRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'simple_demo'-------
    // get current time
    t = simple_demoClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *sd_polygon_abscissa* updates
    if (t >= 0.0 && sd_polygon_abscissa.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sd_polygon_abscissa.tStart = t;  // (not accounting for frame time here)
      sd_polygon_abscissa.frameNStart = frameN;  // exact frame index
      
      sd_polygon_abscissa.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (sd_polygon_abscissa.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      sd_polygon_abscissa.setAutoDraw(false);
    }
    
    // *sd_image_left_seeds* updates
    if (t >= 0.0 && sd_image_left_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sd_image_left_seeds.tStart = t;  // (not accounting for frame time here)
      sd_image_left_seeds.frameNStart = frameN;  // exact frame index
      
      sd_image_left_seeds.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (sd_image_left_seeds.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      sd_image_left_seeds.setAutoDraw(false);
    }
    
    // *sd_image_right_seeds* updates
    if (t >= 0.0 && sd_image_right_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sd_image_right_seeds.tStart = t;  // (not accounting for frame time here)
      sd_image_right_seeds.frameNStart = frameN;  // exact frame index
      
      sd_image_right_seeds.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (sd_image_right_seeds.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      sd_image_right_seeds.setAutoDraw(false);
    }
    
    // *sd_image_pigeon* updates
    if (t >= 0.0 && sd_image_pigeon.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sd_image_pigeon.tStart = t;  // (not accounting for frame time here)
      sd_image_pigeon.frameNStart = frameN;  // exact frame index
      
      sd_image_pigeon.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (sd_image_pigeon.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      sd_image_pigeon.setAutoDraw(false);
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    simple_demoComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function simple_demoRoutineEnd() {
  return async function () {
    //------Ending Routine 'simple_demo'-------
    simple_demoComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    // the Routine "simple_demo" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var status_instruction_setupComponents;
function status_instruction_setupRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'status_instruction_setup'-------
    t = 0;
    status_instruction_setupClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    status_show = true;
    pigeon_show = true;
    seeds_show = true;
    pigeon_flip = 0.1;
    pigeon_true_x = 0;
    pigeon_shown_x = pigeon_true_x;
    box_show = sic_box_show;
    box_width = sic_box_width;
    box_height = sic_box_height;
    box_x = sic_box_x;
    box_y = sic_box_y;
    message_string = sic_message;
    
    // keep track of which components have finished
    status_instruction_setupComponents = [];
    
    status_instruction_setupComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function status_instruction_setupRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'status_instruction_setup'-------
    // get current time
    t = status_instruction_setupClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    status_instruction_setupComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function status_instruction_setupRoutineEnd() {
  return async function () {
    //------Ending Routine 'status_instruction_setup'-------
    status_instruction_setupComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    // the Routine "status_instruction_setup" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var block_name;
var snr_set;
var snr_conditions_file;
var instructions_conditions_file;
var online_trials_count;
var predefined_trials_count;
var start_taskComponents;
function start_taskRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'start_task'-------
    t = 0;
    start_taskClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    trial_start = true;
    trial_number = 1;
    steps_count = 0;
    coins_count = 0;
    last_coins_count = coins_count;
    last_steps_count = steps_count;
    task_type = tc_task_type;
    block_name = tc_block_name;
    steps_max = tc_steps_max;
    steps_taken_to_start_trial = tc_steps_taken_to_start_trial;
    steps_lost_per_error = tc_steps_lost_per_error;
    coins_max = tc_coins_max;
    coins_paid_to_start_trial = tc_coins_paid_to_start_trial;
    coins_lost_per_error = tc_coins_lost_per_error;
    coins_gained_per_correct = tc_coins_gained_per_correct;
    predefined_bound_min_bonus = tc_predefined_bound_min_bonus;
    predefined_bound_max_bonus = tc_predefined_bound_max_bonus;
    snr_set = tc_snr_set;
    if ((snr_set < 10)) {
        snr_conditions_file = (("Resources/PigeonSNRSet0" + snr_set.toString()) + "Conditions.xlsx");
    } else {
        snr_conditions_file = (("Resources/PigeonSNRSet" + snr_set.toString()) + "Conditions.xlsx");
    }
    instructions_conditions_file = (("Resources/Pigeon_" + task_type) + "_InstructionsConditions.xlsx");
    if ((task_type === "online")) {
        online_trials_count = 10000;
        predefined_trials_count = 0;
    } else {
        online_trials_count = 0;
        predefined_trials_count = 10000;
        predefined_bound_x = PREDEFINED_BOUND_X_DEFAULT;
    }
    block_number += 1;
    message_string = (((((("Block " + block_number.toString()) + " of ") + experiment_loop.trialList.length.toString()) + "\n\n") + "Type = ") + task_type);
    pigeon_flip = 0.1;
    pigeon_true_x = 0;
    pigeon_shown_x = pigeon_true_x;
    update_status(true, true, true);
    
    // keep track of which components have finished
    start_taskComponents = [];
    
    start_taskComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function start_taskRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'start_task'-------
    // get current time
    t = start_taskClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    start_taskComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function start_taskRoutineEnd() {
  return async function () {
    //------Ending Routine 'start_task'-------
    start_taskComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    // the Routine "start_task" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var setup_instructions_messageComponents;
function setup_instructions_messageRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'setup_instructions_message'-------
    t = 0;
    setup_instructions_messageClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    message_string = ic_message_string;
    
    // keep track of which components have finished
    setup_instructions_messageComponents = [];
    
    setup_instructions_messageComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function setup_instructions_messageRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'setup_instructions_message'-------
    // get current time
    t = setup_instructions_messageClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    setup_instructions_messageComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function setup_instructions_messageRoutineEnd() {
  return async function () {
    //------Ending Routine 'setup_instructions_message'-------
    setup_instructions_messageComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    // the Routine "setup_instructions_message" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _otr_kbd_response_allKeys;
var correct_ans;
var endRoutine;
var online_task_runComponents;
function online_task_runRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'online_task_run'-------
    t = 0;
    online_task_runClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    otr_kbd_response.keys = undefined;
    otr_kbd_response.rt = undefined;
    _otr_kbd_response_allKeys = [];
    if ((trial_start === true)) {
        otr_kbd_response.clearEvents();
        message_string = "Failed to respond\nPlease respond more quickly";
        instruction_string = CHOICE_STRING;
        pigeon_flip = 0.1;
        pigeon_true_x = ORIGIN_X;
        pigeon_shown_x = pigeon_true_x;
        pigeon_steps = [pigeon_shown_x];
        coins_count += coins_paid_to_start_trial;
        steps_count += steps_taken_to_start_trial;
        if ((util.randint(0, 2) === 1)) {
            correct_ans = "right";
            step_mean_val = Math.abs(step_mean);
        } else {
            correct_ans = "left";
            step_mean_val = ((- 1) * Math.abs(step_mean));
        }
        step_std_val = step_std;
        psychoJS.experiment.addData("step_mean_val", step_mean_val);
        psychoJS.experiment.addData("step_std_val", step_std_val);
        psychoJS.experiment.addData("block_number", block_number);
        psychoJS.experiment.addData("trial_number", trial_number);
        psychoJS.experiment.addData("correct", correct_ans);
        trial_start = false;
    } else {
        steps_count += 1;
        move_pigeon(step_mean_val, step_std_val);
        pigeon_steps.push(pigeon_shown_x);
        update_status(false, false, true);
    }
    if ((steps_count >= steps_max)) {
        steps_count = steps_max;
        psychoJS.experiment.addData("pigeon_steps", pigeon_steps);
        psychoJS.experiment.addData("coins_count", coins_count);
        psychoJS.experiment.addData("steps_count", steps_count);
        psychoJS.experiment.addData("choice", "none");
        set_block_end_string();
        psychoJS.experiment.addData("bonus", bonus_count);
        endRoutine = true;
        online_move_loop.finished = true;
        online_trial_loop.finished = true;
    }
    
    otr_polygon_abscissa.setPos([ORIGIN_X, ABSCISSA_Y]);
    otr_polygon_abscissa.setSize([(EDGE_DISTANCE * 2), LINE_HEIGHT]);
    otr_polygon_ordinate.setPos([ORIGIN_X, ORDINATE_Y]);
    otr_polygon_ordinate.setSize([ORDINATE_WIDTH, ORDINATE_HEIGHT]);
    otr_polygon_petey_midpoint.setPos([pigeon_shown_x, ORDINATE_Y]);
    otr_polygon_petey_midpoint.setSize([(ORDINATE_WIDTH * 0.9), ORDINATE_HEIGHT]);
    otr_image_pigeon.setPos([pigeon_shown_x, ORIGIN_Y]);
    otr_image_pigeon.setSize([pigeon_flip, 0.1]);
    otr_polygon_coins_bar.setPos([COINS_BAR_X, coins_bar_center]);
    otr_polygon_coins_bar.setSize([STATUS_BAR_WIDTH, coins_bar_height]);
    otr_polygon_steps_bar.setPos([STEPS_BAR_X, steps_bar_center]);
    otr_polygon_steps_bar.setSize([STATUS_BAR_WIDTH, steps_bar_height]);
    otr_polygon_status_bar_apex.setPos([ORIGIN_X, STATUS_BAR_APEX_Y]);
    otr_polygon_status_bar_apex.setSize([STATUS_BAR_BRACKET_WIDTH, LINE_HEIGHT]);
    otr_polygon_status_bar_base.setPos([ORIGIN_X, STATUS_BAR_BASE_Y]);
    otr_polygon_status_bar_base.setSize([STATUS_BAR_BRACKET_WIDTH, LINE_HEIGHT]);
    otr_text_coins.setText(coins_count);
    otr_text_coins_properties.setText(coins_string);
    otr_text_steps.setText((steps_max - steps_count));
    otr_text_steps_properties.setText(steps_string);
    otr_text_instruction.setText(instruction_string);
    // keep track of which components have finished
    online_task_runComponents = [];
    online_task_runComponents.push(otr_kbd_response);
    online_task_runComponents.push(otr_polygon_abscissa);
    online_task_runComponents.push(otr_polygon_ordinate);
    online_task_runComponents.push(otr_polygon_petey_midpoint);
    online_task_runComponents.push(otr_image_pigeon);
    online_task_runComponents.push(otr_image_right_seeds);
    online_task_runComponents.push(otr_image_left_seeds);
    online_task_runComponents.push(otr_polygon_coins_bar);
    online_task_runComponents.push(otr_polygon_steps_bar);
    online_task_runComponents.push(otr_polygon_status_bar_apex);
    online_task_runComponents.push(otr_polygon_status_bar_base);
    online_task_runComponents.push(otr_text_coins);
    online_task_runComponents.push(otr_text_coins_label);
    online_task_runComponents.push(otr_text_coins_properties);
    online_task_runComponents.push(otr_text_steps);
    online_task_runComponents.push(otr_text_steps_label);
    online_task_runComponents.push(otr_text_steps_properties);
    online_task_runComponents.push(otr_text_instruction);
    
    online_task_runComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function online_task_runRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'online_task_run'-------
    // get current time
    t = online_task_runClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *otr_kbd_response* updates
    if (t >= 0.0 && otr_kbd_response.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_kbd_response.tStart = t;  // (not accounting for frame time here)
      otr_kbd_response.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { otr_kbd_response.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { otr_kbd_response.start(); }); // start on screen flip
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_kbd_response.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_kbd_response.status = PsychoJS.Status.FINISHED;
  }

    if (otr_kbd_response.status === PsychoJS.Status.STARTED) {
      let theseKeys = otr_kbd_response.getKeys({keyList: ['left', 'right'], waitRelease: false});
      _otr_kbd_response_allKeys = _otr_kbd_response_allKeys.concat(theseKeys);
      if (_otr_kbd_response_allKeys.length > 0) {
        otr_kbd_response.keys = _otr_kbd_response_allKeys[_otr_kbd_response_allKeys.length - 1].name;  // just the last key pressed
        otr_kbd_response.rt = _otr_kbd_response_allKeys[_otr_kbd_response_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *otr_polygon_abscissa* updates
    if (t >= 0 && otr_polygon_abscissa.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_polygon_abscissa.tStart = t;  // (not accounting for frame time here)
      otr_polygon_abscissa.frameNStart = frameN;  // exact frame index
      
      otr_polygon_abscissa.setAutoDraw(true);
    }

    frameRemains = 0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_polygon_abscissa.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_polygon_abscissa.setAutoDraw(false);
    }
    
    // *otr_polygon_ordinate* updates
    if (t >= 0 && otr_polygon_ordinate.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_polygon_ordinate.tStart = t;  // (not accounting for frame time here)
      otr_polygon_ordinate.frameNStart = frameN;  // exact frame index
      
      otr_polygon_ordinate.setAutoDraw(true);
    }

    frameRemains = 0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_polygon_ordinate.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_polygon_ordinate.setAutoDraw(false);
    }
    
    // *otr_polygon_petey_midpoint* updates
    if (t >= 0.0 && otr_polygon_petey_midpoint.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_polygon_petey_midpoint.tStart = t;  // (not accounting for frame time here)
      otr_polygon_petey_midpoint.frameNStart = frameN;  // exact frame index
      
      otr_polygon_petey_midpoint.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_polygon_petey_midpoint.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_polygon_petey_midpoint.setAutoDraw(false);
    }
    
    // *otr_image_pigeon* updates
    if (t >= 0.0 && otr_image_pigeon.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_image_pigeon.tStart = t;  // (not accounting for frame time here)
      otr_image_pigeon.frameNStart = frameN;  // exact frame index
      
      otr_image_pigeon.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_image_pigeon.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_image_pigeon.setAutoDraw(false);
    }
    
    // *otr_image_right_seeds* updates
    if (t >= 0.0 && otr_image_right_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_image_right_seeds.tStart = t;  // (not accounting for frame time here)
      otr_image_right_seeds.frameNStart = frameN;  // exact frame index
      
      otr_image_right_seeds.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_image_right_seeds.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_image_right_seeds.setAutoDraw(false);
    }
    
    // *otr_image_left_seeds* updates
    if (t >= 0.0 && otr_image_left_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_image_left_seeds.tStart = t;  // (not accounting for frame time here)
      otr_image_left_seeds.frameNStart = frameN;  // exact frame index
      
      otr_image_left_seeds.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_image_left_seeds.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_image_left_seeds.setAutoDraw(false);
    }
    
    // *otr_polygon_coins_bar* updates
    if (t >= 0.0 && otr_polygon_coins_bar.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_polygon_coins_bar.tStart = t;  // (not accounting for frame time here)
      otr_polygon_coins_bar.frameNStart = frameN;  // exact frame index
      
      otr_polygon_coins_bar.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_polygon_coins_bar.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_polygon_coins_bar.setAutoDraw(false);
    }
    
    // *otr_polygon_steps_bar* updates
    if (t >= 0.0 && otr_polygon_steps_bar.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_polygon_steps_bar.tStart = t;  // (not accounting for frame time here)
      otr_polygon_steps_bar.frameNStart = frameN;  // exact frame index
      
      otr_polygon_steps_bar.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_polygon_steps_bar.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_polygon_steps_bar.setAutoDraw(false);
    }
    
    // *otr_polygon_status_bar_apex* updates
    if (t >= 0.0 && otr_polygon_status_bar_apex.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_polygon_status_bar_apex.tStart = t;  // (not accounting for frame time here)
      otr_polygon_status_bar_apex.frameNStart = frameN;  // exact frame index
      
      otr_polygon_status_bar_apex.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_polygon_status_bar_apex.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_polygon_status_bar_apex.setAutoDraw(false);
    }
    
    // *otr_polygon_status_bar_base* updates
    if (t >= 0.0 && otr_polygon_status_bar_base.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_polygon_status_bar_base.tStart = t;  // (not accounting for frame time here)
      otr_polygon_status_bar_base.frameNStart = frameN;  // exact frame index
      
      otr_polygon_status_bar_base.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_polygon_status_bar_base.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_polygon_status_bar_base.setAutoDraw(false);
    }
    
    // *otr_text_coins* updates
    if (t >= 0.0 && otr_text_coins.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_text_coins.tStart = t;  // (not accounting for frame time here)
      otr_text_coins.frameNStart = frameN;  // exact frame index
      
      otr_text_coins.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_text_coins.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_text_coins.setAutoDraw(false);
    }
    
    // *otr_text_coins_label* updates
    if (t >= 0.0 && otr_text_coins_label.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_text_coins_label.tStart = t;  // (not accounting for frame time here)
      otr_text_coins_label.frameNStart = frameN;  // exact frame index
      
      otr_text_coins_label.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_text_coins_label.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_text_coins_label.setAutoDraw(false);
    }
    
    // *otr_text_coins_properties* updates
    if (t >= 0.0 && otr_text_coins_properties.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_text_coins_properties.tStart = t;  // (not accounting for frame time here)
      otr_text_coins_properties.frameNStart = frameN;  // exact frame index
      
      otr_text_coins_properties.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_text_coins_properties.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_text_coins_properties.setAutoDraw(false);
    }
    
    // *otr_text_steps* updates
    if (t >= 0.0 && otr_text_steps.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_text_steps.tStart = t;  // (not accounting for frame time here)
      otr_text_steps.frameNStart = frameN;  // exact frame index
      
      otr_text_steps.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_text_steps.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_text_steps.setAutoDraw(false);
    }
    
    // *otr_text_steps_label* updates
    if (t >= 0.0 && otr_text_steps_label.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_text_steps_label.tStart = t;  // (not accounting for frame time here)
      otr_text_steps_label.frameNStart = frameN;  // exact frame index
      
      otr_text_steps_label.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_text_steps_label.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_text_steps_label.setAutoDraw(false);
    }
    
    // *otr_text_steps_properties* updates
    if (t >= 0.0 && otr_text_steps_properties.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_text_steps_properties.tStart = t;  // (not accounting for frame time here)
      otr_text_steps_properties.frameNStart = frameN;  // exact frame index
      
      otr_text_steps_properties.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_text_steps_properties.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_text_steps_properties.setAutoDraw(false);
    }
    
    // *otr_text_instruction* updates
    if (t >= 0.0 && otr_text_instruction.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      otr_text_instruction.tStart = t;  // (not accounting for frame time here)
      otr_text_instruction.frameNStart = frameN;  // exact frame index
      
      otr_text_instruction.setAutoDraw(true);
    }

    frameRemains = 0.0 + UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (otr_text_instruction.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      otr_text_instruction.setAutoDraw(false);
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    online_task_runComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


var _pj;
function online_task_runRoutineEnd() {
  return async function () {
    //------Ending Routine 'online_task_run'-------
    online_task_runComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    otr_kbd_response.stop();
    var _pj;
    function _pj_snippets(container) {
        function in_es6(left, right) {
            if (((right instanceof Array) || ((typeof right) === "string"))) {
                return (right.indexOf(left) > (- 1));
            } else {
                if (((right instanceof Map) || (right instanceof Set) || (right instanceof WeakMap) || (right instanceof WeakSet))) {
                    return right.has(left);
                } else {
                    return (left in right);
                }
            }
        }
        container["in_es6"] = in_es6;
        return container;
    }
    _pj = {};
    _pj_snippets(_pj);
    if (otr_kbd_response.keys) {
        if (_pj.in_es6(correct_ans, otr_kbd_response.keys)) {
            coins_count += coins_gained_per_correct;
            set_trial_end_string(CORRECT_STRING);
            update_status(false, true, false);
        } else {
            coins_count -= coins_lost_per_error;
            steps_count += steps_lost_per_error;
            if ((steps_count >= steps_max)) {
                steps_count = steps_max;
            }
            set_trial_end_string(ERROR_STRING);
            update_status(false, true, true);
        }
        psychoJS.experiment.addData("pigeon_steps", pigeon_steps);
        psychoJS.experiment.addData("coins_count", coins_count);
        psychoJS.experiment.addData("steps_count", steps_count);
        psychoJS.experiment.addData("choice", otr_kbd_response.keys);
        otr_kbd_response.clearEvents();
        trial_number += 1;
        online_move_loop.finished = true;
    }
    
    // the Routine "online_task_run" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _pts_kbd_update_allKeys;
var predefined_task_setupComponents;
function predefined_task_setupRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'predefined_task_setup'-------
    t = 0;
    predefined_task_setupClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    pts_kbd_update.keys = undefined;
    pts_kbd_update.rt = undefined;
    _pts_kbd_update_allKeys = [];
    if (trial_start) {
        trial_start = false;
        pts_kbd_update.clearEvents();
        if ((task_type === "predefined")) {
            instruction_string = PREDEFINED_CHOICE_STRING;
        } else {
            if ((steps_count < min_steps_to_commit)) {
                instruction_string = PREDEFINED_CHOICE_STRING;
            } else {
                instruction_string = PREDEFINED_CHOICE_COMMIT_STRING;
            }
        }
        message_string = PREDEFINED_NO_CHOICE_STRING;
        choice = "max";
        pigeon_flip = 0.1;
        pigeon_true_x = 0;
        pigeon_shown_x = pigeon_true_x;
        pigeon_steps = [pigeon_shown_x];
        predefined_bound_x_previous = predefined_bound_x;
        coins_count += coins_paid_to_start_trial;
        steps_count += steps_taken_to_start_trial;
        step_mean = ssc_step_mean;
        step_std = ssc_step_std;
        if ((util.randint(0, 2) === 1)) {
            correct_ans = "right";
            step_mean_val = Math.abs(step_mean);
        } else {
            correct_ans = "left";
            step_mean_val = ((- 1) * Math.abs(step_mean));
        }
        step_std_val = step_std;
        psychoJS.experiment.addData("step_mean_val", step_mean_val);
        psychoJS.experiment.addData("step_std_val", step_std_val);
        psychoJS.experiment.addData("block_number", block_number);
        psychoJS.experiment.addData("trial_number", trial_number);
        psychoJS.experiment.addData("correct", correct_ans);
        trial_number += 1;
        if ((steps_count >= steps_max)) {
            steps_count = steps_max;
            predefined_setup_loop.finished = true;
            continueRoutine = false;
        }
    }
    
    pts_polygon_abscissa.setPos([ORIGIN_X, ABSCISSA_Y]);
    pts_polygon_abscissa.setSize([(EDGE_DISTANCE * 2), LINE_HEIGHT]);
    pts_polygon_ordinate.setPos([ORIGIN_X, ORDINATE_Y]);
    pts_polygon_ordinate.setSize([ORDINATE_WIDTH, ORDINATE_HEIGHT]);
    pts_polygon_petey_midpoint.setPos([pigeon_shown_x, ORDINATE_Y]);
    pts_polygon_petey_midpoint.setSize([(ORDINATE_WIDTH * 0.9), ORDINATE_HEIGHT]);
    pts_polygon_left_bound.setPos([(- predefined_bound_x), PREDEFINED_BOUND_Y]);
    pts_polygon_left_bound.setSize([PREDEFINED_BOUND_WIDTH, PREDEFINED_BOUND_HEIGHT]);
    pts_polygon_left_bound.setLineColor(new util.Color(PREDEFINED_BOUND_COLOR));
    pts_polygon_right_bound.setPos([predefined_bound_x, PREDEFINED_BOUND_Y]);
    pts_polygon_right_bound.setSize([PREDEFINED_BOUND_WIDTH, PREDEFINED_BOUND_HEIGHT]);
    pts_polygon_stay_left.setPos([(- predefined_bound_x_previous), PREDEFINED_BOUND_MARKER_Y]);
    pts_polygon_stay_left.setSize([PREDEFINED_BOUND_MARKER_SZ, PREDEFINED_BOUND_MARKER_SZ]);
    pts_polygon_stay_right.setPos([predefined_bound_x_previous, PREDEFINED_BOUND_MARKER_Y]);
    pts_polygon_stay_right.setSize([PREDEFINED_BOUND_MARKER_SZ, PREDEFINED_BOUND_MARKER_SZ]);
    pts_image_pigeon.setPos([pigeon_shown_x, ORIGIN_Y]);
    pts_image_pigeon.setSize([pigeon_flip, 0.1]);
    pts_polygon_coins_bar.setPos([COINS_BAR_X, coins_bar_center]);
    pts_polygon_coins_bar.setSize([STATUS_BAR_WIDTH, coins_bar_height]);
    pts_polygon_steps_bar.setPos([STEPS_BAR_X, steps_bar_center]);
    pts_polygon_steps_bar.setSize([STATUS_BAR_WIDTH, steps_bar_height]);
    pts_polygon_status_bar_apex.setPos([ORIGIN_X, STATUS_BAR_APEX_Y]);
    pts_polygon_status_bar_apex.setSize([STATUS_BAR_BRACKET_WIDTH, LINE_HEIGHT]);
    pts_polygon_status_bar_base.setPos([ORIGIN_X, STATUS_BAR_BASE_Y]);
    pts_polygon_status_bar_base.setSize([STATUS_BAR_BRACKET_WIDTH, LINE_HEIGHT]);
    pts_text_coins.setPos([COINS_TEXT_X, COINS_TEXT_Y]);
    pts_text_coins.setText(coins_count);
    pts_text_coins_properties.setText(coins_string);
    pts_text_steps.setText((steps_max - steps_count));
    pts_text_steps_properties.setText(steps_string);
    pts_text_instruction.setText(instruction_string);
    // keep track of which components have finished
    predefined_task_setupComponents = [];
    predefined_task_setupComponents.push(pts_kbd_update);
    predefined_task_setupComponents.push(pts_polygon_abscissa);
    predefined_task_setupComponents.push(pts_polygon_ordinate);
    predefined_task_setupComponents.push(pts_polygon_petey_midpoint);
    predefined_task_setupComponents.push(pts_polygon_left_bound);
    predefined_task_setupComponents.push(pts_polygon_right_bound);
    predefined_task_setupComponents.push(pts_polygon_stay_left);
    predefined_task_setupComponents.push(pts_polygon_stay_right);
    predefined_task_setupComponents.push(pts_image_pigeon);
    predefined_task_setupComponents.push(pts_image_right_seeds);
    predefined_task_setupComponents.push(pts_image_left_seeds);
    predefined_task_setupComponents.push(pts_polygon_coins_bar);
    predefined_task_setupComponents.push(pts_polygon_steps_bar);
    predefined_task_setupComponents.push(pts_polygon_status_bar_apex);
    predefined_task_setupComponents.push(pts_polygon_status_bar_base);
    predefined_task_setupComponents.push(pts_text_coins);
    predefined_task_setupComponents.push(pts_text_coins_label);
    predefined_task_setupComponents.push(pts_text_coins_properties);
    predefined_task_setupComponents.push(pts_text_steps);
    predefined_task_setupComponents.push(pts_text_steps_label);
    predefined_task_setupComponents.push(pts_text_steps_properties);
    predefined_task_setupComponents.push(pts_text_instruction);
    predefined_task_setupComponents.push(pts_text_hint);
    predefined_task_setupComponents.push(pts_text_hint_2);
    
    predefined_task_setupComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function predefined_task_setupRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'predefined_task_setup'-------
    // get current time
    t = predefined_task_setupClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *pts_kbd_update* updates
    if (t >= 0.0 && pts_kbd_update.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_kbd_update.tStart = t;  // (not accounting for frame time here)
      pts_kbd_update.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { pts_kbd_update.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { pts_kbd_update.start(); }); // start on screen flip
    }

    if (pts_kbd_update.status === PsychoJS.Status.STARTED) {
      let theseKeys = pts_kbd_update.getKeys({keyList: ['up', 'down', 'left', 'right', 'space', 'c'], waitRelease: false});
      _pts_kbd_update_allKeys = _pts_kbd_update_allKeys.concat(theseKeys);
      if (_pts_kbd_update_allKeys.length > 0) {
        pts_kbd_update.keys = _pts_kbd_update_allKeys[_pts_kbd_update_allKeys.length - 1].name;  // just the last key pressed
        pts_kbd_update.rt = _pts_kbd_update_allKeys[_pts_kbd_update_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *pts_polygon_abscissa* updates
    if (t >= 0 && pts_polygon_abscissa.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_abscissa.tStart = t;  // (not accounting for frame time here)
      pts_polygon_abscissa.frameNStart = frameN;  // exact frame index
      
      pts_polygon_abscissa.setAutoDraw(true);
    }

    
    // *pts_polygon_ordinate* updates
    if (t >= 0 && pts_polygon_ordinate.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_ordinate.tStart = t;  // (not accounting for frame time here)
      pts_polygon_ordinate.frameNStart = frameN;  // exact frame index
      
      pts_polygon_ordinate.setAutoDraw(true);
    }

    
    // *pts_polygon_petey_midpoint* updates
    if (t >= 0.0 && pts_polygon_petey_midpoint.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_petey_midpoint.tStart = t;  // (not accounting for frame time here)
      pts_polygon_petey_midpoint.frameNStart = frameN;  // exact frame index
      
      pts_polygon_petey_midpoint.setAutoDraw(true);
    }

    
    // *pts_polygon_left_bound* updates
    if (t >= 0.0 && pts_polygon_left_bound.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_left_bound.tStart = t;  // (not accounting for frame time here)
      pts_polygon_left_bound.frameNStart = frameN;  // exact frame index
      
      pts_polygon_left_bound.setAutoDraw(true);
    }

    
    // *pts_polygon_right_bound* updates
    if (t >= 0.0 && pts_polygon_right_bound.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_right_bound.tStart = t;  // (not accounting for frame time here)
      pts_polygon_right_bound.frameNStart = frameN;  // exact frame index
      
      pts_polygon_right_bound.setAutoDraw(true);
    }

    
    // *pts_polygon_stay_left* updates
    if (t >= 0 && pts_polygon_stay_left.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_stay_left.tStart = t;  // (not accounting for frame time here)
      pts_polygon_stay_left.frameNStart = frameN;  // exact frame index
      
      pts_polygon_stay_left.setAutoDraw(true);
    }

    
    // *pts_polygon_stay_right* updates
    if (t >= 0 && pts_polygon_stay_right.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_stay_right.tStart = t;  // (not accounting for frame time here)
      pts_polygon_stay_right.frameNStart = frameN;  // exact frame index
      
      pts_polygon_stay_right.setAutoDraw(true);
    }

    
    // *pts_image_pigeon* updates
    if (t >= 0.0 && pts_image_pigeon.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_image_pigeon.tStart = t;  // (not accounting for frame time here)
      pts_image_pigeon.frameNStart = frameN;  // exact frame index
      
      pts_image_pigeon.setAutoDraw(true);
    }

    
    // *pts_image_right_seeds* updates
    if (t >= 0.0 && pts_image_right_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_image_right_seeds.tStart = t;  // (not accounting for frame time here)
      pts_image_right_seeds.frameNStart = frameN;  // exact frame index
      
      pts_image_right_seeds.setAutoDraw(true);
    }

    
    // *pts_image_left_seeds* updates
    if (t >= 0.0 && pts_image_left_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_image_left_seeds.tStart = t;  // (not accounting for frame time here)
      pts_image_left_seeds.frameNStart = frameN;  // exact frame index
      
      pts_image_left_seeds.setAutoDraw(true);
    }

    
    // *pts_polygon_coins_bar* updates
    if (t >= 0.0 && pts_polygon_coins_bar.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_coins_bar.tStart = t;  // (not accounting for frame time here)
      pts_polygon_coins_bar.frameNStart = frameN;  // exact frame index
      
      pts_polygon_coins_bar.setAutoDraw(true);
    }

    
    // *pts_polygon_steps_bar* updates
    if (t >= 0.0 && pts_polygon_steps_bar.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_steps_bar.tStart = t;  // (not accounting for frame time here)
      pts_polygon_steps_bar.frameNStart = frameN;  // exact frame index
      
      pts_polygon_steps_bar.setAutoDraw(true);
    }

    
    // *pts_polygon_status_bar_apex* updates
    if (t >= 0.0 && pts_polygon_status_bar_apex.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_status_bar_apex.tStart = t;  // (not accounting for frame time here)
      pts_polygon_status_bar_apex.frameNStart = frameN;  // exact frame index
      
      pts_polygon_status_bar_apex.setAutoDraw(true);
    }

    
    // *pts_polygon_status_bar_base* updates
    if (t >= 0.0 && pts_polygon_status_bar_base.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_polygon_status_bar_base.tStart = t;  // (not accounting for frame time here)
      pts_polygon_status_bar_base.frameNStart = frameN;  // exact frame index
      
      pts_polygon_status_bar_base.setAutoDraw(true);
    }

    
    // *pts_text_coins* updates
    if (t >= 0.0 && pts_text_coins.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_text_coins.tStart = t;  // (not accounting for frame time here)
      pts_text_coins.frameNStart = frameN;  // exact frame index
      
      pts_text_coins.setAutoDraw(true);
    }

    
    // *pts_text_coins_label* updates
    if (t >= 0.0 && pts_text_coins_label.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_text_coins_label.tStart = t;  // (not accounting for frame time here)
      pts_text_coins_label.frameNStart = frameN;  // exact frame index
      
      pts_text_coins_label.setAutoDraw(true);
    }

    
    // *pts_text_coins_properties* updates
    if (t >= 0.0 && pts_text_coins_properties.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_text_coins_properties.tStart = t;  // (not accounting for frame time here)
      pts_text_coins_properties.frameNStart = frameN;  // exact frame index
      
      pts_text_coins_properties.setAutoDraw(true);
    }

    
    // *pts_text_steps* updates
    if (t >= 0.0 && pts_text_steps.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_text_steps.tStart = t;  // (not accounting for frame time here)
      pts_text_steps.frameNStart = frameN;  // exact frame index
      
      pts_text_steps.setAutoDraw(true);
    }

    
    // *pts_text_steps_label* updates
    if (t >= 0.0 && pts_text_steps_label.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_text_steps_label.tStart = t;  // (not accounting for frame time here)
      pts_text_steps_label.frameNStart = frameN;  // exact frame index
      
      pts_text_steps_label.setAutoDraw(true);
    }

    
    // *pts_text_steps_properties* updates
    if (t >= 0.0 && pts_text_steps_properties.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_text_steps_properties.tStart = t;  // (not accounting for frame time here)
      pts_text_steps_properties.frameNStart = frameN;  // exact frame index
      
      pts_text_steps_properties.setAutoDraw(true);
    }

    
    // *pts_text_instruction* updates
    if (t >= 0.0 && pts_text_instruction.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_text_instruction.tStart = t;  // (not accounting for frame time here)
      pts_text_instruction.frameNStart = frameN;  // exact frame index
      
      pts_text_instruction.setAutoDraw(true);
    }

    
    // *pts_text_hint* updates
    if (t >= 0.0 && pts_text_hint.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_text_hint.tStart = t;  // (not accounting for frame time here)
      pts_text_hint.frameNStart = frameN;  // exact frame index
      
      pts_text_hint.setAutoDraw(true);
    }

    
    if (pts_text_hint.status === PsychoJS.Status.STARTED){ // only update if being drawn
      pts_text_hint.setPos([ORIGIN_X, MESSAGE_Y], false);
    }
    
    // *pts_text_hint_2* updates
    if (t >= 0.0 && pts_text_hint_2.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      pts_text_hint_2.tStart = t;  // (not accounting for frame time here)
      pts_text_hint_2.frameNStart = frameN;  // exact frame index
      
      pts_text_hint_2.setAutoDraw(true);
    }

    
    if (pts_text_hint_2.status === PsychoJS.Status.STARTED){ // only update if being drawn
      pts_text_hint_2.setPos([ORIGIN_X, 0.2], false);
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    predefined_task_setupComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function predefined_task_setupRoutineEnd() {
  return async function () {
    //------Ending Routine 'predefined_task_setup'-------
    predefined_task_setupComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    pts_kbd_update.stop();
    if ((((task_type === "predefined_commit") && (pts_kbd_update.keys === "c")) && (steps_count > min_steps_to_commit))) {
        psychoJS.experiment.addData("predefined_bound_final", predefined_bound_x);
        psychoJS.experiment.addData("predefined_bound", predefined_bound_x);
        committed_to_bound = true;
        predefined_setup_loop.finished = true;
    } else {
        if ((pts_kbd_update.keys === "space")) {
            psychoJS.experiment.addData("predefined_bound", predefined_bound_x);
            predefined_setup_loop.finished = true;
        } else {
            if (((pts_kbd_update.keys === "left") && (predefined_bound_x <= (EDGE_DISTANCE - PREDEFINED_BOUND_SMALL_DELTA)))) {
                predefined_bound_x += PREDEFINED_BOUND_SMALL_DELTA;
            } else {
                if (((pts_kbd_update.keys === "right") && (predefined_bound_x >= PREDEFINED_BOUND_SMALL_DELTA))) {
                    predefined_bound_x -= PREDEFINED_BOUND_SMALL_DELTA;
                } else {
                    if (((pts_kbd_update.keys === "up") && (predefined_bound_x <= (EDGE_DISTANCE - PREDEFINED_BOUND_LARGE_DELTA)))) {
                        predefined_bound_x += PREDEFINED_BOUND_LARGE_DELTA;
                    } else {
                        if (((pts_kbd_update.keys === "down") && (predefined_bound_x >= PREDEFINED_BOUND_LARGE_DELTA))) {
                            predefined_bound_x -= PREDEFINED_BOUND_LARGE_DELTA;
                        }
                    }
                }
            }
        }
    }
    pts_kbd_update.clearEvents();
    
    // the Routine "predefined_task_setup" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _ptr_kbd_response_allKeys;
var predefined_task_runComponents;
function predefined_task_runRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'predefined_task_run'-------
    t = 0;
    predefined_task_runClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    ptr_kbd_response.keys = undefined;
    ptr_kbd_response.rt = undefined;
    _ptr_kbd_response_allKeys = [];
    if (((committed_to_bound && (steps_count > min_steps_to_commit)) || (steps_count >= steps_max))) {
        if ((steps_count > steps_max)) {
            steps_count = steps_max;
        }
        predefined_move_loop.finished = true;
        predefined_trial_loop.finished = true;
        committed_to_bound = false;
        psychoJS.experiment.addData("steps_count", steps_count);
        steps_count = steps_max;
        update_status(false, false, true);
        set_block_end_string();
        psychoJS.experiment.addData("pigeon_steps", pigeon_steps);
        psychoJS.experiment.addData("coins_count", coins_count);
        psychoJS.experiment.addData("choice", "none");
        psychoJS.experiment.addData("bonus", bonus_count);
        psychoJS.experiment.addData("predefined_bound", predefined_bound_x);
        psychoJS.experiment.addData("predefined_bound_final", predefined_bound_x);
    } else {
        steps_count += 1;
        move_pigeon(step_mean_val, step_std_val);
        pigeon_steps.push(pigeon_shown_x);
        update_status(false, false, true);
        still_going = false;
        if ((pigeon_true_x >= predefined_bound_x)) {
            choice = "right";
        } else {
            if ((pigeon_true_x <= (- predefined_bound_x))) {
                choice = "left";
            } else {
                if ((steps_count >= steps_max)) {
                    choice = "max";
                } else {
                    still_going = true;
                }
            }
        }
        if ((! still_going)) {
            if ((choice === correct_ans)) {
                coins_count += coins_gained_per_correct;
                set_trial_end_string(CORRECT_STRING);
            } else {
                if ((choice !== "max")) {
                    coins_count -= coins_lost_per_error;
                    steps_count += steps_lost_per_error;
                    set_trial_end_string(ERROR_STRING);
                }
            }
            update_status(false, true, false);
            psychoJS.experiment.addData("pigeon_steps", pigeon_steps);
            psychoJS.experiment.addData("coins_count", coins_count);
            psychoJS.experiment.addData("steps_count", steps_count);
            psychoJS.experiment.addData("choice", choice);
            predefined_move_loop.finished = true;
        }
    }
    
    ptr_polygon_abscissa.setPos([ORIGIN_X, ABSCISSA_Y]);
    ptr_polygon_abscissa.setSize([(EDGE_DISTANCE * 2), LINE_HEIGHT]);
    ptr_polygon_ordinate.setPos([ORIGIN_X, ORDINATE_Y]);
    ptr_polygon_ordinate.setSize([ORDINATE_WIDTH, ORDINATE_HEIGHT]);
    ptr_polygon_petey_midpoint.setPos([pigeon_shown_x, ORDINATE_Y]);
    ptr_polygon_petey_midpoint.setSize([(ORDINATE_WIDTH * 0.9), ORDINATE_HEIGHT]);
    ptr_polygon_left_bound.setPos([(- predefined_bound_x), PREDEFINED_BOUND_Y]);
    ptr_polygon_left_bound.setSize([PREDEFINED_BOUND_WIDTH, PREDEFINED_BOUND_HEIGHT]);
    ptr_polygon_left_bound.setLineColor(new util.Color(PREDEFINED_BOUND_COLOR));
    ptr_polygon_right_bound.setPos([predefined_bound_x, PREDEFINED_BOUND_Y]);
    ptr_polygon_right_bound.setSize([PREDEFINED_BOUND_WIDTH, PREDEFINED_BOUND_HEIGHT]);
    ptr_polygon_right_bound.setLineColor(new util.Color(PREDEFINED_BOUND_COLOR));
    ptr_image_pigeon.setPos([pigeon_shown_x, ORIGIN_Y]);
    ptr_image_pigeon.setSize([pigeon_flip, 0.1]);
    ptr_polygon_coins_bar.setPos([COINS_BAR_X, coins_bar_center]);
    ptr_polygon_coins_bar.setSize([STATUS_BAR_WIDTH, coins_bar_height]);
    ptr_polygon_steps_bar.setPos([STEPS_BAR_X, steps_bar_center]);
    ptr_polygon_steps_bar.setSize([STATUS_BAR_WIDTH, steps_bar_height]);
    ptr_polygon_status_bar_apex.setPos([ORIGIN_X, STATUS_BAR_APEX_Y]);
    ptr_polygon_status_bar_apex.setSize([STATUS_BAR_BRACKET_WIDTH, LINE_HEIGHT]);
    ptr_polygon_status_bar_base.setPos([ORIGIN_X, STATUS_BAR_BASE_Y]);
    ptr_polygon_status_bar_base.setSize([STATUS_BAR_BRACKET_WIDTH, LINE_HEIGHT]);
    ptr_text_coins.setText(coins_count);
    ptr_text_coins_properties.setText(coins_string);
    ptr_text_steps.setText((steps_max - steps_count));
    ptr_text_steps_properties.setText(steps_string);
    ptr_text_instruction.setText(PREDEFINED_ABORT_INSTRUCTION);
    // keep track of which components have finished
    predefined_task_runComponents = [];
    predefined_task_runComponents.push(ptr_kbd_response);
    predefined_task_runComponents.push(ptr_polygon_abscissa);
    predefined_task_runComponents.push(ptr_polygon_ordinate);
    predefined_task_runComponents.push(ptr_polygon_petey_midpoint);
    predefined_task_runComponents.push(ptr_polygon_left_bound);
    predefined_task_runComponents.push(ptr_polygon_right_bound);
    predefined_task_runComponents.push(ptr_image_pigeon);
    predefined_task_runComponents.push(ptr_image_right_seeds);
    predefined_task_runComponents.push(ptr_image_left_seeds);
    predefined_task_runComponents.push(ptr_polygon_coins_bar);
    predefined_task_runComponents.push(ptr_polygon_steps_bar);
    predefined_task_runComponents.push(ptr_polygon_status_bar_apex);
    predefined_task_runComponents.push(ptr_polygon_status_bar_base);
    predefined_task_runComponents.push(ptr_text_coins);
    predefined_task_runComponents.push(ptr_text_coins_label);
    predefined_task_runComponents.push(ptr_text_coins_properties);
    predefined_task_runComponents.push(ptr_text_steps);
    predefined_task_runComponents.push(ptr_text_steps_label);
    predefined_task_runComponents.push(ptr_text_steps_properties);
    predefined_task_runComponents.push(ptr_text_instruction);
    
    predefined_task_runComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function predefined_task_runRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'predefined_task_run'-------
    // get current time
    t = predefined_task_runClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *ptr_kbd_response* updates
    if (t >= 0.0 && ptr_kbd_response.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_kbd_response.tStart = t;  // (not accounting for frame time here)
      ptr_kbd_response.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { ptr_kbd_response.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { ptr_kbd_response.start(); }); // start on screen flip
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_kbd_response.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_kbd_response.status = PsychoJS.Status.FINISHED;
  }

    if (ptr_kbd_response.status === PsychoJS.Status.STARTED) {
      let theseKeys = ptr_kbd_response.getKeys({keyList: ['space', 'up', 'down', 'left', 'right'], waitRelease: false});
      _ptr_kbd_response_allKeys = _ptr_kbd_response_allKeys.concat(theseKeys);
      if (_ptr_kbd_response_allKeys.length > 0) {
        ptr_kbd_response.keys = _ptr_kbd_response_allKeys[_ptr_kbd_response_allKeys.length - 1].name;  // just the last key pressed
        ptr_kbd_response.rt = _ptr_kbd_response_allKeys[_ptr_kbd_response_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *ptr_polygon_abscissa* updates
    if (t >= 0 && ptr_polygon_abscissa.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_polygon_abscissa.tStart = t;  // (not accounting for frame time here)
      ptr_polygon_abscissa.frameNStart = frameN;  // exact frame index
      
      ptr_polygon_abscissa.setAutoDraw(true);
    }

    frameRemains = 0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_polygon_abscissa.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_polygon_abscissa.setAutoDraw(false);
    }
    
    // *ptr_polygon_ordinate* updates
    if (t >= 0 && ptr_polygon_ordinate.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_polygon_ordinate.tStart = t;  // (not accounting for frame time here)
      ptr_polygon_ordinate.frameNStart = frameN;  // exact frame index
      
      ptr_polygon_ordinate.setAutoDraw(true);
    }

    frameRemains = 0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_polygon_ordinate.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_polygon_ordinate.setAutoDraw(false);
    }
    
    // *ptr_polygon_petey_midpoint* updates
    if (t >= 0.0 && ptr_polygon_petey_midpoint.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_polygon_petey_midpoint.tStart = t;  // (not accounting for frame time here)
      ptr_polygon_petey_midpoint.frameNStart = frameN;  // exact frame index
      
      ptr_polygon_petey_midpoint.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_polygon_petey_midpoint.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_polygon_petey_midpoint.setAutoDraw(false);
    }
    
    // *ptr_polygon_left_bound* updates
    if (t >= 0.0 && ptr_polygon_left_bound.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_polygon_left_bound.tStart = t;  // (not accounting for frame time here)
      ptr_polygon_left_bound.frameNStart = frameN;  // exact frame index
      
      ptr_polygon_left_bound.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_polygon_left_bound.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_polygon_left_bound.setAutoDraw(false);
    }
    
    // *ptr_polygon_right_bound* updates
    if (t >= 0.0 && ptr_polygon_right_bound.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_polygon_right_bound.tStart = t;  // (not accounting for frame time here)
      ptr_polygon_right_bound.frameNStart = frameN;  // exact frame index
      
      ptr_polygon_right_bound.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_polygon_right_bound.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_polygon_right_bound.setAutoDraw(false);
    }
    
    // *ptr_image_pigeon* updates
    if (t >= 0.0 && ptr_image_pigeon.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_image_pigeon.tStart = t;  // (not accounting for frame time here)
      ptr_image_pigeon.frameNStart = frameN;  // exact frame index
      
      ptr_image_pigeon.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_image_pigeon.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_image_pigeon.setAutoDraw(false);
    }
    
    // *ptr_image_right_seeds* updates
    if (t >= 0.0 && ptr_image_right_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_image_right_seeds.tStart = t;  // (not accounting for frame time here)
      ptr_image_right_seeds.frameNStart = frameN;  // exact frame index
      
      ptr_image_right_seeds.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_image_right_seeds.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_image_right_seeds.setAutoDraw(false);
    }
    
    // *ptr_image_left_seeds* updates
    if (t >= 0.0 && ptr_image_left_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_image_left_seeds.tStart = t;  // (not accounting for frame time here)
      ptr_image_left_seeds.frameNStart = frameN;  // exact frame index
      
      ptr_image_left_seeds.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_image_left_seeds.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_image_left_seeds.setAutoDraw(false);
    }
    
    // *ptr_polygon_coins_bar* updates
    if (t >= 0.0 && ptr_polygon_coins_bar.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_polygon_coins_bar.tStart = t;  // (not accounting for frame time here)
      ptr_polygon_coins_bar.frameNStart = frameN;  // exact frame index
      
      ptr_polygon_coins_bar.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_polygon_coins_bar.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_polygon_coins_bar.setAutoDraw(false);
    }
    
    // *ptr_polygon_steps_bar* updates
    if (t >= 0.0 && ptr_polygon_steps_bar.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_polygon_steps_bar.tStart = t;  // (not accounting for frame time here)
      ptr_polygon_steps_bar.frameNStart = frameN;  // exact frame index
      
      ptr_polygon_steps_bar.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_polygon_steps_bar.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_polygon_steps_bar.setAutoDraw(false);
    }
    
    // *ptr_polygon_status_bar_apex* updates
    if (t >= 0.0 && ptr_polygon_status_bar_apex.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_polygon_status_bar_apex.tStart = t;  // (not accounting for frame time here)
      ptr_polygon_status_bar_apex.frameNStart = frameN;  // exact frame index
      
      ptr_polygon_status_bar_apex.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_polygon_status_bar_apex.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_polygon_status_bar_apex.setAutoDraw(false);
    }
    
    // *ptr_polygon_status_bar_base* updates
    if (t >= 0.0 && ptr_polygon_status_bar_base.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_polygon_status_bar_base.tStart = t;  // (not accounting for frame time here)
      ptr_polygon_status_bar_base.frameNStart = frameN;  // exact frame index
      
      ptr_polygon_status_bar_base.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_polygon_status_bar_base.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_polygon_status_bar_base.setAutoDraw(false);
    }
    
    // *ptr_text_coins* updates
    if (t >= 0.0 && ptr_text_coins.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_text_coins.tStart = t;  // (not accounting for frame time here)
      ptr_text_coins.frameNStart = frameN;  // exact frame index
      
      ptr_text_coins.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_text_coins.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_text_coins.setAutoDraw(false);
    }
    
    // *ptr_text_coins_label* updates
    if (t >= 0.0 && ptr_text_coins_label.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_text_coins_label.tStart = t;  // (not accounting for frame time here)
      ptr_text_coins_label.frameNStart = frameN;  // exact frame index
      
      ptr_text_coins_label.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_text_coins_label.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_text_coins_label.setAutoDraw(false);
    }
    
    // *ptr_text_coins_properties* updates
    if (t >= 0.0 && ptr_text_coins_properties.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_text_coins_properties.tStart = t;  // (not accounting for frame time here)
      ptr_text_coins_properties.frameNStart = frameN;  // exact frame index
      
      ptr_text_coins_properties.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_text_coins_properties.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_text_coins_properties.setAutoDraw(false);
    }
    
    // *ptr_text_steps* updates
    if (t >= 0.0 && ptr_text_steps.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_text_steps.tStart = t;  // (not accounting for frame time here)
      ptr_text_steps.frameNStart = frameN;  // exact frame index
      
      ptr_text_steps.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_text_steps.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_text_steps.setAutoDraw(false);
    }
    
    // *ptr_text_steps_label* updates
    if (t >= 0.0 && ptr_text_steps_label.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_text_steps_label.tStart = t;  // (not accounting for frame time here)
      ptr_text_steps_label.frameNStart = frameN;  // exact frame index
      
      ptr_text_steps_label.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_text_steps_label.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_text_steps_label.setAutoDraw(false);
    }
    
    // *ptr_text_steps_properties* updates
    if (t >= 0.0 && ptr_text_steps_properties.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_text_steps_properties.tStart = t;  // (not accounting for frame time here)
      ptr_text_steps_properties.frameNStart = frameN;  // exact frame index
      
      ptr_text_steps_properties.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_text_steps_properties.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_text_steps_properties.setAutoDraw(false);
    }
    
    // *ptr_text_instruction* updates
    if (t >= 0.0 && ptr_text_instruction.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      ptr_text_instruction.tStart = t;  // (not accounting for frame time here)
      ptr_text_instruction.frameNStart = frameN;  // exact frame index
      
      ptr_text_instruction.setAutoDraw(true);
    }

    frameRemains = 0.0 + PREDEFINED_UPDATE_INTERVAL - psychoJS.window.monitorFramePeriod * 0.75;  // most of one frame period left
    if (ptr_text_instruction.status === PsychoJS.Status.STARTED && t >= frameRemains) {
      ptr_text_instruction.setAutoDraw(false);
    }
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    predefined_task_runComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function predefined_task_runRoutineEnd() {
  return async function () {
    //------Ending Routine 'predefined_task_run'-------
    predefined_task_runComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    ptr_kbd_response.stop();
    if (ptr_kbd_response.keys) {
        ptr_kbd_response.clearEvents();
        choice = "none";
        set_trial_end_string(PREDEFINED_ABORT_STRING);
        predefined_move_loop.finished = true;
    }
    psychoJS.experiment.addData("pigeon_steps", pigeon_steps);
    psychoJS.experiment.addData("coins_count", coins_count);
    psychoJS.experiment.addData("steps_count", steps_count);
    psychoJS.experiment.addData("choice", choice);
    
    // the Routine "predefined_task_run" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var final_messageComponents;
function final_messageRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'final_message'-------
    t = 0;
    final_messageClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    message_string = (("Finished, fantastic job!\nYou earned " + bonus_count.toString()) + " bonuses.");
    
    // keep track of which components have finished
    final_messageComponents = [];
    
    final_messageComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function final_messageRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'final_message'-------
    // get current time
    t = final_messageClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    final_messageComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function final_messageRoutineEnd() {
  return async function () {
    //------Ending Routine 'final_message'-------
    final_messageComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    // the Routine "final_message" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


var _sfm_kbd_continue_allKeys;
var show_final_messageComponents;
function show_final_messageRoutineBegin(snapshot) {
  return async function () {
    TrialHandler.fromSnapshot(snapshot); // ensure that .thisN vals are up to date
    
    //------Prepare to start Routine 'show_final_message'-------
    t = 0;
    show_final_messageClock.reset(); // clock
    frameN = -1;
    continueRoutine = true; // until we're told otherwise
    // update component parameters for each repeat
    sfm_kbd_continue.keys = undefined;
    sfm_kbd_continue.rt = undefined;
    _sfm_kbd_continue_allKeys = [];
    sm_kbd_continue.clearEvents();
    
    sfm_text_message.setText(message_string);
    sfm_text_steps_properties.setText(steps_string);
    sfm_text_steps.setText((steps_max - steps_count));
    sfm_text_coins_properties.setText(coins_string);
    sfm_text_coins.setText(coins_count);
    sfm_polygon_status_bar_base.setPos([ORIGIN_X, STATUS_BAR_BASE_Y]);
    sfm_polygon_status_bar_base.setSize([STATUS_BAR_BRACKET_WIDTH, LINE_HEIGHT]);
    sfm_polygon_status_bar_apex.setPos([ORIGIN_X, STATUS_BAR_APEX_Y]);
    sfm_polygon_status_bar_apex.setSize([STATUS_BAR_BRACKET_WIDTH, LINE_HEIGHT]);
    sfm_polygon_steps_bar.setPos([STEPS_BAR_X, steps_bar_center]);
    sfm_polygon_steps_bar.setSize([STATUS_BAR_WIDTH, steps_bar_height]);
    sfm_polygon_coin_bar.setPos([COINS_BAR_X, coins_bar_center]);
    sfm_polygon_coin_bar.setSize([STATUS_BAR_WIDTH, coins_bar_height]);
    sfm_image_petey.setPos([pigeon_shown_x, ORIGIN_Y]);
    sfm_image_petey.setSize([pigeon_flip, 0.1]);
    sfm_polygon_petey_midpoint.setPos([pigeon_shown_x, ORDINATE_Y]);
    sfm_polygon_petey_midpoint.setSize([(ORDINATE_WIDTH * 0.9), ORDINATE_HEIGHT]);
    sfm_polygon_ordinate.setPos([ORIGIN_X, ORDINATE_Y]);
    sfm_polygon_ordinate.setSize([ORDINATE_WIDTH, ORDINATE_HEIGHT]);
    sfm_polygon_abscissa.setPos([ORIGIN_X, ABSCISSA_Y]);
    sfm_polygon_abscissa.setSize([ABSCISSA_WIDTH, LINE_HEIGHT]);
    sfm_polygon_box.setPos([box_x, box_y]);
    sfm_polygon_box.setSize([box_width, box_height]);
    sfm_text_instruction.setText(CONTINUE_STRING);
    // keep track of which components have finished
    show_final_messageComponents = [];
    show_final_messageComponents.push(sfm_kbd_continue);
    show_final_messageComponents.push(sfm_text_message);
    show_final_messageComponents.push(sfm_text_steps_properties);
    show_final_messageComponents.push(sfm_text_steps_label);
    show_final_messageComponents.push(sfm_text_steps);
    show_final_messageComponents.push(sfm_text_coins_properties);
    show_final_messageComponents.push(sfm_text_coins_label);
    show_final_messageComponents.push(sfm_text_coins);
    show_final_messageComponents.push(sfm_polygon_status_bar_base);
    show_final_messageComponents.push(sfm_polygon_status_bar_apex);
    show_final_messageComponents.push(sfm_polygon_steps_bar);
    show_final_messageComponents.push(sfm_polygon_coin_bar);
    show_final_messageComponents.push(sfm_image_petey);
    show_final_messageComponents.push(sfm_image_left_seeds);
    show_final_messageComponents.push(sfm_image_right_seeds);
    show_final_messageComponents.push(sfm_polygon_petey_midpoint);
    show_final_messageComponents.push(sfm_polygon_ordinate);
    show_final_messageComponents.push(sfm_polygon_abscissa);
    show_final_messageComponents.push(sfm_polygon_box);
    show_final_messageComponents.push(sfm_text_instruction);
    
    show_final_messageComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent)
        thisComponent.status = PsychoJS.Status.NOT_STARTED;
       });
    return Scheduler.Event.NEXT;
  }
}


function show_final_messageRoutineEachFrame() {
  return async function () {
    //------Loop for each frame of Routine 'show_final_message'-------
    // get current time
    t = show_final_messageClock.getTime();
    frameN = frameN + 1;// number of completed frames (so 0 is the first frame)
    // update/draw components on each frame
    
    // *sfm_kbd_continue* updates
    if (t >= 0.5 && sfm_kbd_continue.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_kbd_continue.tStart = t;  // (not accounting for frame time here)
      sfm_kbd_continue.frameNStart = frameN;  // exact frame index
      
      // keyboard checking is just starting
      psychoJS.window.callOnFlip(function() { sfm_kbd_continue.clock.reset(); });  // t=0 on next screen flip
      psychoJS.window.callOnFlip(function() { sfm_kbd_continue.start(); }); // start on screen flip
      psychoJS.window.callOnFlip(function() { sfm_kbd_continue.clearEvents(); });
    }

    if (sfm_kbd_continue.status === PsychoJS.Status.STARTED) {
      let theseKeys = sfm_kbd_continue.getKeys({keyList: ['space', 'up', 'down', 'left', 'right'], waitRelease: false});
      _sfm_kbd_continue_allKeys = _sfm_kbd_continue_allKeys.concat(theseKeys);
      if (_sfm_kbd_continue_allKeys.length > 0) {
        sfm_kbd_continue.keys = _sfm_kbd_continue_allKeys[_sfm_kbd_continue_allKeys.length - 1].name;  // just the last key pressed
        sfm_kbd_continue.rt = _sfm_kbd_continue_allKeys[_sfm_kbd_continue_allKeys.length - 1].rt;
        // a response ends the routine
        continueRoutine = false;
      }
    }
    
    
    // *sfm_text_message* updates
    if ((message_show) && sfm_text_message.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_text_message.tStart = t;  // (not accounting for frame time here)
      sfm_text_message.frameNStart = frameN;  // exact frame index
      
      sfm_text_message.setAutoDraw(true);
    }

    
    // *sfm_text_steps_properties* updates
    if ((status_show) && sfm_text_steps_properties.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_text_steps_properties.tStart = t;  // (not accounting for frame time here)
      sfm_text_steps_properties.frameNStart = frameN;  // exact frame index
      
      sfm_text_steps_properties.setAutoDraw(true);
    }

    
    // *sfm_text_steps_label* updates
    if ((status_show) && sfm_text_steps_label.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_text_steps_label.tStart = t;  // (not accounting for frame time here)
      sfm_text_steps_label.frameNStart = frameN;  // exact frame index
      
      sfm_text_steps_label.setAutoDraw(true);
    }

    
    // *sfm_text_steps* updates
    if ((status_show) && sfm_text_steps.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_text_steps.tStart = t;  // (not accounting for frame time here)
      sfm_text_steps.frameNStart = frameN;  // exact frame index
      
      sfm_text_steps.setAutoDraw(true);
    }

    
    // *sfm_text_coins_properties* updates
    if ((status_show) && sfm_text_coins_properties.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_text_coins_properties.tStart = t;  // (not accounting for frame time here)
      sfm_text_coins_properties.frameNStart = frameN;  // exact frame index
      
      sfm_text_coins_properties.setAutoDraw(true);
    }

    
    // *sfm_text_coins_label* updates
    if ((status_show) && sfm_text_coins_label.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_text_coins_label.tStart = t;  // (not accounting for frame time here)
      sfm_text_coins_label.frameNStart = frameN;  // exact frame index
      
      sfm_text_coins_label.setAutoDraw(true);
    }

    
    // *sfm_text_coins* updates
    if ((status_show) && sfm_text_coins.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_text_coins.tStart = t;  // (not accounting for frame time here)
      sfm_text_coins.frameNStart = frameN;  // exact frame index
      
      sfm_text_coins.setAutoDraw(true);
    }

    
    // *sfm_polygon_status_bar_base* updates
    if ((status_show) && sfm_polygon_status_bar_base.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_polygon_status_bar_base.tStart = t;  // (not accounting for frame time here)
      sfm_polygon_status_bar_base.frameNStart = frameN;  // exact frame index
      
      sfm_polygon_status_bar_base.setAutoDraw(true);
    }

    
    // *sfm_polygon_status_bar_apex* updates
    if ((status_show) && sfm_polygon_status_bar_apex.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_polygon_status_bar_apex.tStart = t;  // (not accounting for frame time here)
      sfm_polygon_status_bar_apex.frameNStart = frameN;  // exact frame index
      
      sfm_polygon_status_bar_apex.setAutoDraw(true);
    }

    
    // *sfm_polygon_steps_bar* updates
    if ((status_show) && sfm_polygon_steps_bar.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_polygon_steps_bar.tStart = t;  // (not accounting for frame time here)
      sfm_polygon_steps_bar.frameNStart = frameN;  // exact frame index
      
      sfm_polygon_steps_bar.setAutoDraw(true);
    }

    
    // *sfm_polygon_coin_bar* updates
    if ((status_show) && sfm_polygon_coin_bar.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_polygon_coin_bar.tStart = t;  // (not accounting for frame time here)
      sfm_polygon_coin_bar.frameNStart = frameN;  // exact frame index
      
      sfm_polygon_coin_bar.setAutoDraw(true);
    }

    
    // *sfm_image_petey* updates
    if ((pigeon_show) && sfm_image_petey.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_image_petey.tStart = t;  // (not accounting for frame time here)
      sfm_image_petey.frameNStart = frameN;  // exact frame index
      
      sfm_image_petey.setAutoDraw(true);
    }

    
    // *sfm_image_left_seeds* updates
    if ((seeds_show) && sfm_image_left_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_image_left_seeds.tStart = t;  // (not accounting for frame time here)
      sfm_image_left_seeds.frameNStart = frameN;  // exact frame index
      
      sfm_image_left_seeds.setAutoDraw(true);
    }

    
    // *sfm_image_right_seeds* updates
    if ((seeds_show) && sfm_image_right_seeds.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_image_right_seeds.tStart = t;  // (not accounting for frame time here)
      sfm_image_right_seeds.frameNStart = frameN;  // exact frame index
      
      sfm_image_right_seeds.setAutoDraw(true);
    }

    
    // *sfm_polygon_petey_midpoint* updates
    if ((status_show) && sfm_polygon_petey_midpoint.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_polygon_petey_midpoint.tStart = t;  // (not accounting for frame time here)
      sfm_polygon_petey_midpoint.frameNStart = frameN;  // exact frame index
      
      sfm_polygon_petey_midpoint.setAutoDraw(true);
    }

    
    // *sfm_polygon_ordinate* updates
    if ((status_show) && sfm_polygon_ordinate.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_polygon_ordinate.tStart = t;  // (not accounting for frame time here)
      sfm_polygon_ordinate.frameNStart = frameN;  // exact frame index
      
      sfm_polygon_ordinate.setAutoDraw(true);
    }

    
    // *sfm_polygon_abscissa* updates
    if ((pigeon_show) && sfm_polygon_abscissa.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_polygon_abscissa.tStart = t;  // (not accounting for frame time here)
      sfm_polygon_abscissa.frameNStart = frameN;  // exact frame index
      
      sfm_polygon_abscissa.setAutoDraw(true);
    }

    
    // *sfm_polygon_box* updates
    if ((box_show) && sfm_polygon_box.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_polygon_box.tStart = t;  // (not accounting for frame time here)
      sfm_polygon_box.frameNStart = frameN;  // exact frame index
      
      sfm_polygon_box.setAutoDraw(true);
    }

    
    // *sfm_text_instruction* updates
    if ((instruction_show) && sfm_text_instruction.status === PsychoJS.Status.NOT_STARTED) {
      // keep track of start time/frame for later
      sfm_text_instruction.tStart = t;  // (not accounting for frame time here)
      sfm_text_instruction.frameNStart = frameN;  // exact frame index
      
      sfm_text_instruction.setAutoDraw(true);
    }

    // check for quit (typically the Esc key)
    if (psychoJS.experiment.experimentEnded || psychoJS.eventManager.getKeys({keyList:['escape']}).length > 0) {
      return quitPsychoJS('The [Escape] key was pressed. Goodbye!', false);
    }
    
    // check if the Routine should terminate
    if (!continueRoutine) {  // a component has requested a forced-end of Routine
      return Scheduler.Event.NEXT;
    }
    
    continueRoutine = false;  // reverts to True if at least one component still running
    show_final_messageComponents.forEach( function(thisComponent) {
      if ('status' in thisComponent && thisComponent.status !== PsychoJS.Status.FINISHED) {
        continueRoutine = true;
      }
    });
    
    // refresh the screen if continuing
    if (continueRoutine) {
      return Scheduler.Event.FLIP_REPEAT;
    } else {
      return Scheduler.Event.NEXT;
    }
  };
}


function show_final_messageRoutineEnd() {
  return async function () {
    //------Ending Routine 'show_final_message'-------
    show_final_messageComponents.forEach( function(thisComponent) {
      if (typeof thisComponent.setAutoDraw === 'function') {
        thisComponent.setAutoDraw(false);
      }
    });
    sfm_kbd_continue.stop();
    trial_start = true;
    
    // the Routine "show_final_message" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset();
    
    return Scheduler.Event.NEXT;
  };
}


function endLoopIteration(scheduler, snapshot) {
  // ------Prepare for next entry------
  return async function () {
    if (typeof snapshot !== 'undefined') {
      // ------Check if user ended loop early------
      if (snapshot.finished) {
        // Check for and save orphaned data
        if (psychoJS.experiment.isEntryEmpty()) {
          psychoJS.experiment.nextEntry(snapshot);
        }
        scheduler.stop();
      } else {
        const thisTrial = snapshot.getCurrentTrial();
        if (typeof thisTrial === 'undefined' || !('isTrials' in thisTrial) || thisTrial.isTrials) {
          psychoJS.experiment.nextEntry(snapshot);
        }
      }
    return Scheduler.Event.NEXT;
    }
  };
}


function importConditions(currentLoop) {
  return async function () {
    psychoJS.importAttributes(currentLoop.getCurrentTrial());
    return Scheduler.Event.NEXT;
    };
}


async function quitPsychoJS(message, isCompleted) {
  // Check for and save orphaned data
  if (psychoJS.experiment.isEntryEmpty()) {
    psychoJS.experiment.nextEntry();
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  psychoJS.window.close();
  psychoJS.quit({message: message, isCompleted: isCompleted});
  
  return Scheduler.Event.QUIT;
}
