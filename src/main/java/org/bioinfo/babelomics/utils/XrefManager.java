package org.bioinfo.babelomics.utils;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.sun.jersey.api.client.Client;
import com.sun.jersey.api.client.ClientResponse;
import com.sun.jersey.api.client.WebResource;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;

/**
 * Created by ralonso on 12/26/14.
 */
public class XrefManager {

    private List<String> list;
    private String species;
    private String dbname;
    private int requests;
    private Map<String, List<String>> listXref;

    private String uri;
    private String babelomicsHome = System.getenv("BABELOMICS_HOME");


//    public XrefManager(List<String> list, String species) {
//        this(list, species, "ensembl_transcript");
//
//    }


    public XrefManager(List<String> list, String species, String dbname) {
        this.list = list;
        this.species = species;
        this.dbname = dbname;
        this.listXref = new HashMap<String, List<String>>();
        this.process();

    }

    private void process() {
        Properties prop = new Properties();
        InputStream input = null;

        try {
            String file = this.babelomicsHome + "/conf/webservices.conf";
            input = new FileInputStream(file);
            prop.load(input);
            this.uri = prop.getProperty("HOST");
            String auxSpecies = species;

            if (species.equalsIgnoreCase("hsa"))
                auxSpecies = "hsapiens";

            this.uri += auxSpecies + "/feature/id/";

            this.requests = Integer.parseInt(prop.getProperty("REQUESTS"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Get map with this structure: [featureId:anotationsList] *
     */


    public Map<String, List<String>> getXrefs() {
        int numberIds = 0;
        StringBuilder batch = new StringBuilder();
        for (String id : list) {
            if (numberIds == requests) {
                this.fillXref(batch);
                numberIds = 0;
                batch = new StringBuilder();
            }
            batch.append(id);
            batch.append(",");
            numberIds++;
        }
        this.fillXref(batch);
        return listXref;
    }

    public void fillXref(StringBuilder batch) {
        Client client = Client.create();

        int lastIndexOfComma = batch.lastIndexOf(",");
        batch.deleteCharAt(lastIndexOfComma);

        String uri = this.uri + batch.toString() + "/xref?dbname=" + this.dbname;
        System.out.println("uri = " + uri);
//                WebResource webResource = client.resource("https://www.ebi.ac.uk/cellbase/webservices/rest/v3/hsapiens/feature/id/" + batch.toString() + "/xref?dbname=ensembl_transcript");
        WebResource webResource = client.resource(uri);
        ClientResponse response = webResource.get(ClientResponse.class);
        String resp = response.getEntity(String.class);
        ObjectMapper mapper = new ObjectMapper();
        JsonNode actualObj = null;
        try {
            actualObj = mapper.readTree(resp).get("response");
        } catch (IOException e) {
            e.printStackTrace();
        }
        for (JsonNode jnode : actualObj) {
            if (listXref.containsKey(jnode.get("id")))
                continue;
            List<String> ids = new ArrayList<String>();
            listXref.put(jnode.get("id").asText(), ids);
            Iterator<JsonNode> it = jnode.get("result").iterator();
            while (it.hasNext()) {
                JsonNode node = it.next();
                ids.add(node.get("id").asText());
            }
        }
    }
}
