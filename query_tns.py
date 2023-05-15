# Perform a search on the Transient Name Server
# https://www.wis-tns.org/
def tns_search(search_obj):
    search_url = url_tns_api + "/search"
    headers = {
        "User-Agent": 'tns_marker{"tns_id": "'
        + str(TNS_BOT_ID)
        + '", "type": "bot", "name": "'
        + TNS_BOT_NAME
        + '"}'
    }
    json_file = OrderedDict(search_obj)
    search_data = {"api_key": TNS_API_KEY, "data": json.dumps(json_file)}
    response = requests.post(search_url, headers=headers, data=search_data)
    return response


# Download the info of the object of interest from the TNS
def get_tns(get_obj):
    get_url = url_tns_api + "/object"
    headers = {
        "User-Agent": 'tns_marker{"tns_id": "'
        + str(TNS_BOT_ID)
        + '", "type": "bot", "name": "'
        + TNS_BOT_NAME
        + '"}'
    }
    json_file = OrderedDict(get_obj)
    get_data = {"api_key": TNS_API_KEY, "data": json.dumps(json_file)}
    response = requests.post(get_url, headers=headers, data=get_data)
    return response
