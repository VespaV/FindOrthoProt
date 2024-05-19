from django.urls import path
from . import views
from . import tree_visualisation

urlpatterns = [
    path('', views.analysis_seq_view, name='sequence_analysis'),
    path('analysis/', views.analysis_success_view, name='analysis_success'),
    path('analysis/tree_visualisation/', tree_visualisation.handle_request, name='handle_request')
]