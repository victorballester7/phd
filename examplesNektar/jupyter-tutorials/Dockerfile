FROM dmoxey/nektar-workbook:20211014_1

RUN pip install --no-cache-dir nbgitpuller

USER root
COPY ./entrypoint.sh /
RUN chmod +x /entrypoint.sh

USER $NB_USER
ENTRYPOINT ["/entrypoint.sh"]


