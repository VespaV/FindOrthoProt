
Чтобы запустить проект нужно ввести команду python3 manage.py runserver сайт будет доступен по адресу: http://127.0.0.1:8000
Необходимые пакеты:
pip3 install daphne
pip3 install channels
pip3 install django-bootstrap-v5
pip3 install Bio
python3 -m pip3 install Django>
Чтобы запустить websocket: daphne -p 8081 FindOrthoProt.asgi:application

python3 manage.py runserver 8081
daphne -p 8081 FindOrthoProt.asgi:application

daphne -p 8080 FindOrthoProt.asgi:application


MONGO DB
to start mogo: brew services start mongodb-community
path: mongodb://localhost:27017/

to start celary
celery -A FindOrthoProt worker -l info    
celery -A FindOrthoProt beat -l info     




