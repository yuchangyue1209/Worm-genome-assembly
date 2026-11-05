#!/usr/bin/env bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
MODULE_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
CONFIG=${CONFIG:-"${MODULE_DIR}/config.sh"}

if [[ ! -s "${CONFIG}" ]]; then
    echo "ERROR: missing ${CONFIG}; copy config.example.sh to config.sh" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "${CONFIG}"

require_file() {
    [[ -s "$1" ]] || { echo "ERROR: missing or empty $1" >&2; exit 1; }
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "ERROR: command not found: $1" >&2
        exit 1
    }
}
