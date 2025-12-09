#!/bin/bash

# ApisTox Deployment Verification Script
# Tests all endpoints after Railway deployment completes

set -e

BASE_URL="https://web-production-f014a.up.railway.app"
TIMEOUT=10

echo "======================================"
echo "ApisTox Deployment Verification"
echo "======================================"
echo ""
echo "Base URL: $BASE_URL"
echo "Timeout: ${TIMEOUT}s per request"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test 1: Health Check
echo -e "${YELLOW}[1/5] Testing Health Endpoint${NC}"
echo "GET $BASE_URL/health"
HEALTH_RESPONSE=$(curl -s -w "\n%{http_code}" --max-time $TIMEOUT "$BASE_URL/health" || echo "")
HEALTH_CODE=$(echo "$HEALTH_RESPONSE" | tail -1)
HEALTH_BODY=$(echo "$HEALTH_RESPONSE" | head -n -1)

if [ "$HEALTH_CODE" = "200" ]; then
    echo -e "${GREEN}✓ Health check passed${NC}"
    echo "$HEALTH_BODY" | jq '.' 2>/dev/null || echo "$HEALTH_BODY"
else
    echo -e "${RED}✗ Health check failed (HTTP $HEALTH_CODE)${NC}"
    echo "$HEALTH_BODY"
fi
echo ""

# Test 2: Root Endpoint
echo -e "${YELLOW}[2/5] Testing Root Endpoint${NC}"
echo "GET $BASE_URL/"
ROOT_RESPONSE=$(curl -s -w "\n%{http_code}" --max-time $TIMEOUT "$BASE_URL/" || echo "")
ROOT_CODE=$(echo "$ROOT_RESPONSE" | tail -1)
ROOT_BODY=$(echo "$ROOT_RESPONSE" | head -n -1)

if [ "$ROOT_CODE" = "200" ]; then
    echo -e "${GREEN}✓ Root endpoint passed${NC}"
    echo "$ROOT_BODY" | jq '.' 2>/dev/null || echo "$ROOT_BODY"
else
    echo -e "${RED}✗ Root endpoint failed (HTTP $ROOT_CODE)${NC}"
    echo "$ROOT_BODY"
fi
echo ""

# Test 3: Test Endpoint (Built-in Sample)
echo -e "${YELLOW}[3/5] Testing Built-in Sample Prediction${NC}"
echo "GET $BASE_URL/test"
TEST_RESPONSE=$(curl -s -w "\n%{http_code}" --max-time $TIMEOUT "$BASE_URL/test" || echo "")
TEST_CODE=$(echo "$TEST_RESPONSE" | tail -1)
TEST_BODY=$(echo "$TEST_RESPONSE" | head -n -1)

if [ "$TEST_CODE" = "200" ]; then
    echo -e "${GREEN}✓ Test endpoint passed${NC}"
    echo "$TEST_BODY" | jq '.prediction' 2>/dev/null || echo "$TEST_BODY"
else
    echo -e "${RED}✗ Test endpoint failed (HTTP $TEST_CODE)${NC}"
    echo "$TEST_BODY"
fi
echo ""

# Test 4: Prediction Endpoint
echo -e "${YELLOW}[4/5] Testing Prediction Endpoint${NC}"
echo "POST $BASE_URL/predict"
PREDICT_DATA='{
  "name": "Imidacloprid",
  "smiles": "C1=CN=C(N1)NC(=O)NCCl",
  "category": "Insecticide",
  "molecular_weight": 255.66,
  "logp": 0.57,
  "exposure_route": "Contact"
}'
echo "Payload: $PREDICT_DATA"
PREDICT_RESPONSE=$(curl -s -w "\n%{http_code}" --max-time $TIMEOUT \
  -X POST "$BASE_URL/predict" \
  -H "Content-Type: application/json" \
  -d "$PREDICT_DATA" || echo "")
PREDICT_CODE=$(echo "$PREDICT_RESPONSE" | tail -1)
PREDICT_BODY=$(echo "$PREDICT_RESPONSE" | head -n -1)

if [ "$PREDICT_CODE" = "200" ]; then
    echo -e "${GREEN}✓ Prediction endpoint passed${NC}"
    echo "$PREDICT_BODY" | jq '.' 2>/dev/null || echo "$PREDICT_BODY"
else
    echo -e "${RED}✗ Prediction endpoint failed (HTTP $PREDICT_CODE)${NC}"
    echo "$PREDICT_BODY"
fi
echo ""

# Test 5: API Documentation
echo -e "${YELLOW}[5/5] Testing API Documentation${NC}"
echo "GET $BASE_URL/docs"
DOCS_RESPONSE=$(curl -s -w "\n%{http_code}" --max-time $TIMEOUT "$BASE_URL/docs" || echo "")
DOCS_CODE=$(echo "$DOCS_RESPONSE" | tail -1)

if [ "$DOCS_CODE" = "200" ]; then
    echo -e "${GREEN}✓ API docs accessible${NC}"
    echo "Visit: $BASE_URL/docs"
else
    echo -e "${RED}✗ API docs failed (HTTP $DOCS_CODE)${NC}"
fi
echo ""

# Summary
echo "======================================"
echo "Summary"
echo "======================================"
PASS_COUNT=0
[ "$HEALTH_CODE" = "200" ] && PASS_COUNT=$((PASS_COUNT+1))
[ "$ROOT_CODE" = "200" ] && PASS_COUNT=$((PASS_COUNT+1))
[ "$TEST_CODE" = "200" ] && PASS_COUNT=$((PASS_COUNT+1))
[ "$PREDICT_CODE" = "200" ] && PASS_COUNT=$((PASS_COUNT+1))
[ "$DOCS_CODE" = "200" ] && PASS_COUNT=$((PASS_COUNT+1))

echo "Tests Passed: $PASS_COUNT/5"
echo ""

if [ "$PASS_COUNT" -eq 5 ]; then
    echo -e "${GREEN}🎉 All tests passed! Deployment is successful.${NC}"
    exit 0
elif [ "$PASS_COUNT" -eq 0 ]; then
    echo -e "${RED}❌ All tests failed. Deployment may still be building or has issues.${NC}"
    echo ""
    echo "Troubleshooting:"
    echo "1. Check Railway dashboard: https://railway.app/dashboard"
    echo "2. Wait 10 minutes from git push time"
    echo "3. Review deployment logs for errors"
    exit 1
else
    echo -e "${YELLOW}⚠️  Some tests passed, some failed. Check logs for details.${NC}"
    exit 1
fi
