#!/bin/bash

# Initialize variables with default values
event_id=""
magnitude=""
longitude=""
latitude=""
depth=""

server_endpoint="risk.crs.ogs.it:8080"

# Function to display usage information
usage() {
    echo "Usage: $0 <event_id> -m <magnitude> -lon <longitude> -lat <latitude> [-d <depth>]"
    exit 1
}

# Check if the event ID is provided
if [ "$#" -lt 5 ]; then
    usage
fi

# Get the first parameter as the event ID
event_id="$1"
shift

# Parse command line options
while [ "$#" -gt 0 ]; do
    case "$1" in
        -m|--magnitude)
            shift
            magnitude="$1"
            ;;
        -lon|--longitude)
            shift
            longitude="$1"
            ;;
        -lat|--latitude)
            shift
            latitude="$1"
            ;;
        -d|--depth)
            shift
            depth="$1"
            ;;
        *)
            usage
            ;;
    esac
    shift
done

# Check if all required parameters are provided
if [ -z "$magnitude" ] || [ -z "$longitude" ] || [ -z "$latitude" ]; then
    usage
fi

# Data to send to the server
data="{\"event_id\": \"$event_id\", \"magnitude\": $magnitude, \"longitude\": $longitude, \"latitude\": $latitude"

# Include depth in data if provided
if [ -n "$depth" ]; then
    data="$data, \"depth\": $depth"
fi

data="$data}"

# Make the HTTP POST request using curl
#curl -X POST -H "Content-Type: application/json" -d "$data" "$server_endpoint"

echo "$data"
