import logging
import os

def setup_logging():
    log_file = "loggest.txt"
    if not os.path.exists(log_file):
        open(log_file, "w").close()  # Crea el archivo vacío si no existe

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )