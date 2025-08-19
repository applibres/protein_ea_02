# Imagen base con Python
FROM rosettacommons/rosetta

# Establece el directorio de trabajo
WORKDIR /app

# Copia todo
COPY . .

RUN apt-get update && apt-get install -y nano
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

CMD ["bash"]
