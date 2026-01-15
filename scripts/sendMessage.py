import http.client
from urllib.parse import urlencode
import credentialsTelegram as ct

TOKEN = ct.TOKEN
CHAT_ID = ct.CHAT_ID

def send_telegram_message(message, details=""):
    print("Sending message to Telegram...")
    host = "api.telegram.org"
    endpoint = f"/bot{TOKEN}/sendMessage"
    conn = http.client.HTTPSConnection(host)
    headers = {"Content-Type": "application/x-www-form-urlencoded"}

    try:
        # Define the API endpoint and parameters
        params = {
            "chat_id": CHAT_ID,
            "text": message + details,
            "parse_mode": "HTML",
            "disable_notification": True
        }
        # Encode the parameters
        body = urlencode(params)
        conn.request("POST", endpoint, body, headers)

        # Get the response
        response = conn.getresponse()
        if response.status == 200:
            print("Message sent successfully!")
        else:
            print(f"Failed to send message: {response.status} {response.reason}")
            print(f"Response: {response.read().decode()}")
            print("Try only without details")
            
            # Retry without details
            params = {
                "chat_id": CHAT_ID,
                "text": message,
                "parse_mode": "HTML",
                "disable_notification": True
            }

            # Encode the parameters
            body = urlencode(params)
            conn.request("POST", endpoint, body, headers)
            response = conn.getresponse()
            if response.status == 200:
                print("Message sent successfully without details!")
            else:
                print(f"Failed to send message without details: {response.status} {response.reason}")
                print(f"Response: {response.read().decode()}")
    except Exception as e:
        print(f"An error occurred while sending the message: {e}")

