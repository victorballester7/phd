import os
import socket
import sendMessage as sm
import sys

def notify_program_start(job_id):
    current_directory = os.getcwd()
    hostname = socket.gethostname()

    message = (
        "✅ <b>Program started</b> 🚀\n\n"
        f"📂 <b>Directory</b>: <code>{current_directory}</code>\n"
        f"💻 <b>Hostname</b>: <code>{hostname}</code>\n"
        f"🆔 <b>Job ID</b>: <code>{job_id}</code>"
    )

    try:
        sm.send_telegram_message(message)
    except Exception as e:
        print(f"Failed to send notification: {e}")

if __name__ == "__main__":
    job_id = sys.argv[1]
    notify_program_start(job_id)

