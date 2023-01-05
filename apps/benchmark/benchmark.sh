#!/bin/bash
#make sure the background processes are killed if the script is interrupted:
trap "kill 0" EXIT 
#print help function:
help () {
    echo "This is the benchmarking tool for CECILIA. It can simulate latency, jitter and bandwidth restrictions."
    echo
    echo "positional arguments:"
    echo "1: number of operations (i.e. vector length for vectorised operations)"
    echo "2: vector length (for operations on vectors)"
    echo "3 & 4: matrix size (x, y) (gram matrix has size x*x) (for operations on matrices, not vectorised vector operations)"
    echo "5: window size (used in CNN)"
    echo "6: kernel size (used in CNN)"
    echo "7: number of kernel (used in CNN)"
    echo "8: number of cycles (how often to run each function for the benchmark)"
    echo "(9+: which functions to run (e.g. MUL, MUX; MAXPOOL for MAX on windows); without this, all functions are run)"
    echo
    echo "options (with defaults included):"
    echo "  -h,   --help            show this help message and exit"
    echo "  -m=d, --mode=d          select debug (d) or release (r) build"
    echo "  -l=0, --latency=0       add a latency of 0 ms to each request"
    echo "  -j=0, --jitter=0        randomly vary the latency by +/-0 ms (cannot be larger than latency)"
    echo "  -b=∞, --bandwidth=∞     cap the bandwidth at ∞ KB/s"
}
# parse arguments:
JITTER=0
MODE=d
while [[ $# -gt 0 ]]; do
    i=$1
    case $i in
        -h|--help)
            help
            exit 0
            ;;
        -l=*|--latency=*)
            LATENCY="${i#*=}"
            shift
            ;;
        -j=*|--jitter=*)
            JITTER="${i#*=}"
            shift
            ;;
        -b=*|--bandwidth=*)
            BANDWIDTH="${i#*=}"
            shift
            ;;
        -m=*|--mode=*)
            MODE="${i#*=}"
            shift
            ;;
        -*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            POSITIONAL_ARGS+=("$1") # save positional arg
            shift
            ;;
    esac
done
# make sure that enough positional arguments were provided:
if [[ ${#POSITIONAL_ARGS[@]} -lt 8 ]]
then
  echo "Positional arguments were not provided. Call this tool with -h to see which ones are required."
  exit 1
fi
# make positional arguments into string:
ARGS=${POSITIONAL_ARGS[*]}
# select appropriate executable:
if [[ "$MODE" = d ]]
then
  EXECUTABLE=../../cmake-build-debug/benchmark
elif [[ "$MODE" = r ]]
then
  EXECUTABLE=../../cmake-build-release/benchmark
else
  echo "Unknown mode \"$MODE\". Proceeding with default (debug)."
fi
# make sure that latency is larger than jitter:
if [[ "$JITTER" -gt 0 ]]
then
  if ! [[ "$LATENCY" ]] || [[ "$LATENCY" -lt "$JITTER" ]]
  then
    echo "Jitter cannot be higher than Latency."
    exit 1
  fi
fi
PORT_P1P2=6381
PORT_P2P1=6381
PORT_HELPERP=6379
PORT_PHELPER=6379
# set up toxiproxy if it is used:
if [[ "$LATENCY" ]] || [[ "$BANDWIDTH" ]]
then
  PORT_P2P1=26381
  PORT_PHELPER=26379
  # find toxiproxy executables:
  TOXI_SERVER_ARRAY=(toxiproxy-server*)
  TOXI_SERVER=${TOXI_SERVER_ARRAY[0]}
  TOXI_CLI_ARRAY=(toxiproxy-cli*)
  TOXI_CLI=${TOXI_CLI_ARRAY[0]}
  # initialise toxiproxy:
  "./$TOXI_SERVER" &>toxi.log &
  sleep 1
  "./$TOXI_CLI" create \
  -l localhost:"$PORT_PHELPER" \
  -u localhost:"$PORT_HELPERP" \
  helper_p1 &>toxi.log
  "./$TOXI_CLI" create \
  -l localhost:"$((PORT_PHELPER + 1))" \
  -u localhost:"$((PORT_HELPERP + 1))" \
  helper_p2 &>toxi.log
  "./$TOXI_CLI" create \
  -l localhost:"$PORT_P2P1" \
  -u localhost:"$PORT_P1P2" \
  p1_p2 &>toxi.log
  if [[ "$LATENCY" ]]
  then
    "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --downstream helper_p1 &>toxi.log
    "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --upstream helper_p1 &>toxi.log
    "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --downstream helper_p2 &>toxi.log
    "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --upstream helper_p2 &>toxi.log
    "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --downstream p1_p2 &>toxi.log
    "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --upstream p1_p2 &>toxi.log
  fi
  if [[ "$BANDWIDTH" ]]
  then
    "./$TOXI_CLI" toxic add -t bandwith -a rate="$BANDWIDTH" --downstream helper_p1 &>toxi.log
    "./$TOXI_CLI" toxic add -t bandwith -a rate="$BANDWIDTH" --upstream helper_p1 &>toxi.log
    "./$TOXI_CLI" toxic add -t bandwith -a rate="$BANDWIDTH" --downstream helper_p2 &>toxi.log
    "./$TOXI_CLI" toxic add -t bandwith -a rate="$BANDWIDTH" --upstream helper_p2 &>toxi.log
    "./$TOXI_CLI" toxic add -t bandwith -a rate="$BANDWIDTH" --downstream p1_p2 &>toxi.log
    "./$TOXI_CLI" toxic add -t bandwith -a rate="$BANDWIDTH" --upstream p1_p2 &>toxi.log
  fi
fi
# run benchmarks:
$EXECUTABLE HELPER $PORT_HELPERP 127.0.0.1 0 0 $ARGS &
sleep 2
$EXECUTABLE P1 $PORT_PHELPER 127.0.0.1 $PORT_P1P2 127.0.0.1 $ARGS &
sleep 2
$EXECUTABLE P2 $PORT_PHELPER 127.0.0.1 $PORT_P2P1 127.0.0.1 $ARGS

