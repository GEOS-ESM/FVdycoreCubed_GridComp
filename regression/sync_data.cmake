# Retrieve temporary AWS credentials from the SI-team credentials API
execute_process(
  COMMAND curl -s -H "x-api-key: $ENV{AWS_API_KEY_GMAO_SITEAM_S3}"
          https://llvsm4u7ij.execute-api.us-east-1.amazonaws.com/credentials
  OUTPUT_VARIABLE CREDS_JSON
  RESULT_VARIABLE CURL_RESULT
)
if(CURL_RESULT)
  message(FATAL_ERROR "Failed to retrieve AWS credentials (curl exit code: ${CURL_RESULT})")
endif()

string(JSON AWS_ACCESS_KEY_ID     GET "${CREDS_JSON}" AccessKeyId)
string(JSON AWS_SECRET_ACCESS_KEY GET "${CREDS_JSON}" SecretAccessKey)
string(JSON AWS_SESSION_TOKEN     GET "${CREDS_JSON}" SessionToken)

# Sync regression test data from S3
execute_process(
  COMMAND ${CMAKE_COMMAND} -E env
    AWS_EC2_METADATA_DISABLED=true
    AWS_ACCESS_KEY_ID=${AWS_ACCESS_KEY_ID}
    AWS_SECRET_ACCESS_KEY=${AWS_SECRET_ACCESS_KEY}
    AWS_SESSION_TOKEN=${AWS_SESSION_TOKEN}
    ${AWS_CLI} s3 sync --only-show-errors "${S3_URI}" "${LOCAL_DIR}"
  RESULT_VARIABLE AWS_RESULT
)
if(AWS_RESULT)
  message(FATAL_ERROR "aws s3 sync failed (exit code: ${AWS_RESULT})")
endif()
