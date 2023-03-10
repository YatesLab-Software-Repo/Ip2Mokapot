FROM python:3.10-bullseye

WORKDIR /usr/src/app

COPY requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD streamlit run ./home.py --server.maxUploadSize 2000