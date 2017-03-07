// This file contains functionality for interfacing with the
// compiled GABE javascript

function init()
{
  var params = {
    "nx" : 128,
    "ny" : 128,
    "nz" : 128,
    "steps" : 100,
    "plot_interval" : 20, // milliseconds
    "plot_data" : [{
       z: [[0,0],[0,0]],
       type: 'surface'
    }],
    "layout" : {
      title: 'Simulation data',
      autosize: true,
      width: 500,
      height: 500,
      margin: {
        l: 65,
        r: 50,
        b: 65,
        t: 90,
      }
    }, // layout
  }; // params

  initializePlot(params);
  simInit(params);
  plotResults(params);
} // init()

async function runSteps(steps)
{
  for(s = 0; s <= steps; s++)
  {
    console.log("Running step " + s + "/" + steps);
    await sleep(50); // overhead to slow things down...
    runStep();
  }
}

function getField(nx, ny, nz)
{
  copyout_fld = Module.cwrap(
    'copyout_fld', 'number', ['number', 'number', 'number']
  );

  var fld = new Float32Array(new Array(nx*ny*nz).fill(0));

  var nDataBytes = fld.length * fld.BYTES_PER_ELEMENT;
  var dataPtr = Module._malloc(nDataBytes);
  var dataHeap = new Uint8Array(Module.HEAPU8.buffer, dataPtr, nDataBytes);
  dataHeap.set(new Uint8Array(fld.buffer));

  copyout_fld(0, dataHeap.byteOffset, fld.length);
  var result = new Float32Array(dataHeap.buffer, dataHeap.byteOffset, fld.length);

  Module._free(dataHeap.byteOffset);
  return result;
}

function simInit(params)
{
  sim_init = Module.cwrap(
    'sim_init', 'number',
    ['number', 'number', 'number']
  );
  var msg = sim_init(params.nx, params.ny, params.nz);
  if(msg == 7) {
    console.log("Simulation initialized!");
  } else {
    console.log("Error initializing sim!");
  }
  return msg;
}

function runStep()
{
   sim_step = Module.cwrap(
    'sim_step', 'number',
    ['number', 'number', 'number']
  );
  var msg = sim_step();
  if(msg == 7) {
    console.log("Simulation step run!");
  } else {
    console.log("Error running sim step!");
  }
  return msg;
}

function initializePlot(params)
{
  params.graphDiv = $("#graphDiv")[0];
  Plotly.newPlot('graphDiv', params.plot_data, params.layout);
}

// asynchronous function, calls itself to continuously plot results
async function plotResults(params)
{
  console.log("Plotting results...");
  var data = getField(params.nx, params.ny, params.nz);

  // get a 2-d slice of data
  var z_data = [];
  for(i=0; i<params.nx; i++) {
    z_data[i] = [];
    for(j=0; j<params.ny; j++) {
      z_data[i][j] = data[i*params.ny*params.nz + j*params.nz + Math.floor(params.nz/2.)];
    }
  }

  // do this instead of making a new plot each time:
  // https://plot.ly/javascript/plotlyjs-function-reference/#plotly-redraw
  params.graphDiv.data[0].z = z_data;
  Plotly.redraw(graphDiv);

  await sleep(params.plot_interval);
  plotResults(params);
}

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}
