from django.urls import path


from . import views

urlpatterns = [
    path('', views.homologous_sequence_view, name='index'),
]