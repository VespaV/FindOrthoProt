from django.urls import path
from . import views

urlpatterns = [
    path('', views.analysis_seq_view, name='sequence_analysis'),
    path('analysis/', views.analysis_success_view, name='analysis_success'),
    path('meme/', views.meme_view, name='meme')
]