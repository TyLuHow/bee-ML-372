import { GoogleGenAI, Type, SchemaType } from "@google/genai";
import { ChemicalData, PredictionResult } from "../types";

const ai = new GoogleGenAI({ apiKey: process.env.API_KEY });

export const analyzeChemicalToxicity = async (data: ChemicalData): Promise<PredictionResult> => {
  try {
    const prompt = `
      Act as an expert ecotoxicologist specializing in Apis mellifera (honey bees).
      Analyze the following chemical compound properties to predict toxicity risk:
      
      Compound Name/Type: ${data.name || "Unknown generic pesticide"}
      Molecular Weight: ${data.mw} g/mol
      LogP (Lipophilicity): ${data.logP}
      Exposure Route: ${data.exposure}
      Category: ${data.category}

      Provide a risk assessment based on general chemical structure-activity relationships (QSAR) principles for bees.
      
      Return the response in strict JSON format.
    `;

    const response = await ai.models.generateContent({
      model: "gemini-2.5-flash",
      contents: prompt,
      config: {
        responseMimeType: "application/json",
        responseSchema: {
          type: Type.OBJECT,
          properties: {
            toxicity: {
              type: Type.STRING,
              enum: ["Toxic", "Safe", "Uncertain"],
              description: "The predicted toxicity classification.",
            },
            confidence: {
              type: Type.NUMBER,
              description: "Confidence score between 0 and 100.",
            },
            explanation: {
              type: Type.STRING,
              description: "A 2-sentence scientific explanation of why.",
            },
            recommendation: {
              type: Type.STRING,
              description: "One actionable recommendation for researchers.",
            },
          },
          required: ["toxicity", "confidence", "explanation", "recommendation"],
        },
      },
    });

    if (response.text) {
      return JSON.parse(response.text) as PredictionResult;
    }

    throw new Error("No response text");
  } catch (error) {
    console.error("Gemini Analysis Error:", error);
    // Fallback mock response if API fails or key is missing
    return {
      toxicity: "Uncertain",
      confidence: 0,
      explanation: "Unable to connect to AI analysis engine. Please verify API key.",
      recommendation: "Proceed with standard in vivo bioassays.",
    };
  }
};
