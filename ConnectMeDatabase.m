function dbc = ConnectMeDatabase()
% Put connection settings here:

connSettings = struct();
connSettings.hostname = 'localhost';
connSettings.dbname = 'GODaily';
connSettings.username = 'root';
connSettings.password = 'ben1234';
dbc = SQL_opendatabase(connSettings);

end
