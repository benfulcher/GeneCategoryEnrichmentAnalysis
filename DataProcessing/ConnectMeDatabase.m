function dbc = ConnectMeDatabase()
% ConnectMetDatabase    Form a connection to a mySQL database
%-------------------------------------------------------------------------------
% PLACE YOUR MYSQL CONNECTION SETTINGS HERE:
%-------------------------------------------------------------------------------

connSettings = struct();
connSettings.hostname = 'localhost';
connSettings.dbname = 'GODaily';
connSettings.username = 'root';
connSettings.password = 'ben1234';
dbc = SQL_opendatabase(connSettings);

end
