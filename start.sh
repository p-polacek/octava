#!/bin/bash

if [ ! -e data/.noupdate ]; then
  # Pull the latest code from GitHub
  echo "Pulling the latest code from GitHub"

  if git pull origin main; then
    echo "Update done"
  else
    echo "Update failed"
    touch updatefailed
  fi
fi

# Set the FLASK_APP environment variable
export FLASK_APP=flask_app.py

# Start the Flask app
echo "Starting"
flask run --host=0.0.0.0 --port=5101