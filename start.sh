#!/bin/bash

# Pull the latest code from GitHub
git reset --hard
echo "Pulling the latest code from GitHub"
git pull || touch .noupdate

roudgpi || touch .noupdate

chmod +x start.sh

# Set the FLASK_APP environment variable
export FLASK_APP=flask_app.py

# Start the Flask app
echo "Starting"
flask run --host=0.0.0.0 --port=5101