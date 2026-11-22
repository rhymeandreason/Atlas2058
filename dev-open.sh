#!/usr/bin/env bash
# Serve the static site at http://localhost:8765 (detached).
cd "$(dirname "$0")" || exit 1
PORT=8765
if lsof -ti tcp:$PORT >/dev/null 2>&1; then
  echo "already running on localhost:$PORT"
  exit 0
fi
nohup python3 -m http.server $PORT --bind 127.0.0.1 >/tmp/atlas2058-dev.log 2>&1 &
echo $! > .dev.pid
echo "serving http://localhost:$PORT"
