"""
URL configuration for FindOrthoProt project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import include, path
from FindOrthoProt import check_format
from django.apps import apps
from django.core.management import call_command


urlpatterns = [
    path('admin/', admin.site.urls),
    path('sequence_analysis/', include('sequence_analysis.urls')),
    path('check_format/', check_format.check_format_view, name='check_format'),
    path('', include('homologous_sequence.urls')),
]


if __name__ == '__main__':
    if apps.is_installed('django_celery_beat'):
        call_command('celery', '-A', 'FindOrthoProt', 'worker', '-l', 'info')

