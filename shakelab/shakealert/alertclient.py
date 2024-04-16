import argparse
import requests
import json

def send_data(event_id, magnitude, longitude, latitude, depth=None,
              server_endpoint="risk.crs.ogs.it:8080"):

    # Construct the data payload
    data = {
        "event_id": event_id,
        "magnitude": magnitude,
        "longitude": longitude,
        "latitude": latitude
    }

    # Include depth in data if provided
    if depth is not None:
        data["depth"] = depth

    # Convert data to JSON format
    json_data = json.dumps(data)

    # Make the HTTP POST request using requests
    #response = requests.post(server_endpoint,
    #                         headers={"Content-Type": "application/json"},
    #                         data=json_data)

    # Print the response
    #print(response.text)
    print(json_data)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Send data to a server")

    parser.add_argument("event_id", help="Event ID")

    parser.add_argument("-m", "--magnitude", type=float,
                        required=True, help="Magnitude")

    parser.add_argument("-lon", "--longitude", type=float,
                        required=True, help="Longitude")

    parser.add_argument("-lat", "--latitude", type=float,
                        required=True, help="Latitude")

    parser.add_argument("-d", "--depth", type=float,
                        required=False, help="Depth (optional)")
    
    args = parser.parse_args()

    send_data(args.event_id,
              args.magnitude,
              args.longitude,
              args.latitude,
              args.depth)
