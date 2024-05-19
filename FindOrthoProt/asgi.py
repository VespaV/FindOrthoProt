import os
from django.core.asgi import get_asgi_application
from channels.routing import ProtocolTypeRouter, URLRouter
from channels.auth import AuthMiddlewareStack
import FindOrthoProt.routing  # Импорт вашего модуля с роутингом

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FindOrthoProt.settings')

application = ProtocolTypeRouter(
    {
        "http": get_asgi_application(),
        "websocket": AuthMiddlewareStack(
            URLRouter(
                FindOrthoProt.routing.websocket_urlpatterns
            )
        ),
    }
)