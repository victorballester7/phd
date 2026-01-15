import sys
import sendMessage as sm
import os
import socket

MAX_CHAR_TELEGRAM_MSG = 4096


def get_last_lines(file, num_lines=40, scale=0.8):
    last_lines = []
    with open(file, "r") as f:
        lines = f.readlines()
        num_characters = MAX_CHAR_TELEGRAM_MSG + 1
        while num_characters > MAX_CHAR_TELEGRAM_MSG:
            last_lines = lines[-int(num_lines * scale) :]
            num_characters = sum(len(line) for line in last_lines)
            scale *= scale
    return last_lines


def get_first_lines(file, num_lines=40, scale=0.8):
    MAX_CHAR_TELEGRAM_MSG = 4096
    first_lines = []
    with open(file, "r") as f:
        lines = f.readlines()
        num_characters = MAX_CHAR_TELEGRAM_MSG + 1
        while num_characters > MAX_CHAR_TELEGRAM_MSG:
            first_lines = lines[: int(num_lines * scale)]
            num_characters = sum(len(line) for line in first_lines)
            scale *= scale
    return first_lines


def check_and_notify(job_id):
    output_file = "output.txt"
    log_file = "log.txt"

    current_directory = os.getcwd()
    hostname = socket.gethostname()

    try:
        log_lines = get_last_lines(log_file)
        output_lines = get_last_lines(output_file)
        timeout = any(("timeout" in line.lower()) for line in log_lines)
        finished = any("Elapsed time Avg" in line for line in output_lines)
        print(f"Timeout: {timeout}, Finished: {finished}")
        if timeout:
            message = (
                "⚠️ <b>Program timed out</b> ⏳\n\n"
                    f"📂 <b>Directory</b>: <code>{current_directory}</code>\n"
                    f"💻 <b>Hostname</b>: <code>{hostname}</code>\n"
                    f"🆔 <b>Job ID</b>: <code>{job_id}</code>"
            )
            try:
                details = (
                    "\n\n<b>Output so far:</b>\n<blockquote expandable>...\n"
                        + "".join(output_lines)
                        + "</blockquote>"
                )
            except Exception:
                details = f"\n\n No <code>{output_file}</code> file found."
        elif not log_lines or finished:
            message = (
                "✅ <b>Program finished successfully</b> 🏁\n\n"
                    f"📂 <b>Directory</b>: <code>{current_directory}</code>\n"
                    f"💻 <b>Hostname</b>: <code>{hostname}</code>\n"
                    f"🆔 <b>Job ID</b>: <code>{job_id}</code>"
            )
            try:
                details = (
                    "\n\n<b>Output:</b>\n<blockquote expandable>...\n"
                        + "".join(output_lines)
                        + "</blockquote>"
                )
            except Exception:
                details = f"\n\n No <code>{output_file}</code> file found."
        else:
            message = (
                "❌ <b>Program encountered an error</b> ⛔\n\n"
                f"📂 <b>Directory</b>: <code>{current_directory}</code>\n"
                f"💻 <b>Hostname</b>: <code>{hostname}</code>\n"
                f"🆔 <b>Job ID</b>: <code>{job_id}</code>"
            )
            details = (
                "\n\n<b>Error Log:</b>\n<blockquote expandable>"
                + "".join(log_lines)
                + "</blockquote>"
            )

        try:
            sm.send_telegram_message(message, details)
        except Exception as e:
            with open("log.txt", "a") as log_file:
                log_file.write(f"Failed to send notification: {e}")
    except Exception as e:
        error_message = (
            "⚠️ <b>An error occurred while checking log files.</b>\n\n"
            f"📂 <b>Directory</b>: <code>{current_directory}</code>\n"
            f"💻 <b>Hostname</b>: <code>{hostname}</code>\n"
            f"🆔 <b>Job ID</b>: <code>{job_id}</code>\n\n"
            f"<b>Error:</b> <code>{e}</code>"
        )
        try:
            sm.send_telegram_message(error_message)
        except Exception as e:
            with open("log.txt", "a") as log_file:
                log_file.write(f"Failed to send notification: {e}")


if __name__ == "__main__":
    job_id = sys.argv[1]
    check_and_notify(job_id)
