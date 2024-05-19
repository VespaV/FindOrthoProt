
from django.contrib import admin
from django.urls import include, path
from FindOrthoProt import check_format
from sequence_analysis import tree_visualisation

urlpatterns = [
    path('admin/', admin.site.urls),
    path('sequence_analysis/', include('sequence_analysis.urls')),
    path('check_format/', check_format.check_format_view, name='check_format'),
    path('', include('homologous_sequence.urls')),
]

