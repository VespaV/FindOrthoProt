from channels.routing import ProtocolTypeRouter, URLRouter
from django.urls import re_path
from .consumers import MyConsumer

websocket_urlpatterns = [
    re_path(r'ws/some_path/$', MyConsumer.as_asgi()),
]

application = ProtocolTypeRouter(
    {
        "websocket": URLRouter(
            websocket_urlpatterns
        ),
    }
)