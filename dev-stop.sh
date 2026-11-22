#!/usr/bin/env bash
# Stop the static server started by dev-open.sh.
cd "$(dirname "$0")" || exit 1
PORT=8765
if [ -f .dev.pid ]; then
  kill "$(cat .dev.pid)" 2>/dev/null
  rm -f .dev.pid
fi
# Catch a server started some other way on the same port.
lsof -ti tcp:$PORT -sTCP:LISTEN | xargs kill 2>/dev/null
exit 0
