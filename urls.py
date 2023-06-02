from django.urls import path
from django.conf.urls import url


from .views import (
index,
get_structures,
documentation,
download_table,
download_session,
)


urlpatterns = [
    path('', index, name='index'),
    url(r'^structure$', get_structures, name='get_structures'),
    url(r'^structure/(?P<pdb>[a-zA-Z0-9.,]{4})$', get_structures, name='get_structures'),
    url(r'^download_table/', download_table, name='download_table'),
    url(r'^download_session/(?P<pdb>^[a-zA-Z0-9]{4})$', download_session, name='download_session'),
    url(r'^download_session/(?P<filename>.*?.pdb)', download_session, name='download_session'),
    url(r'^documentation', documentation, name='documentation'),
]