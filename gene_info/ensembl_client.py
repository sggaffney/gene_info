# https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client

import requests
import time

# class RequestFailedException(Exception):
#     pass

class LookupFailedException(Exception):
    pass

class EnsemblRestClient(object):
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=dict(), params=dict(),
                            data_dict=dict()):

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        # if params:
        #     endpoint += '?' + urllib.urlencode(params)

        out_data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        if not data_dict:
            r = requests.get(self.server + endpoint,
                             headers=hdrs,
                             params=params)
        else:
            r = requests.post(self.server + endpoint,
                              headers=hdrs,
                              params=params,
                              json=data_dict)

        self.req_count += 1

        try:
            r.raise_for_status()
            out_data = r.json()

        except requests.exceptions.HTTPError as e:
            # check if we are being rate limited by the server
            if e.response.status_code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    print("Retrying after {}".format(retry))
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                raise LookupFailedException(
                    'Request failed for {0}: Status code: {1.status_code} '
                    'Reason: {1.reason}\n'.format(endpoint, e.response))

        return out_data
