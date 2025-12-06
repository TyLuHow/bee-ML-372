<objective>
Remove the Gemini API key dependency from the ApisTox application and replace it with mock data for demonstration purposes. This was a styling mockup, not a functional AI integration, so the app should work without requiring any API keys.
</objective>

<context>
The ApisTox bee toxicity prediction platform currently requires a GEMINI_API_KEY environment variable to function. This is causing deployment issues with a white screen error: "An API Key must be set when running in a browser".

The original purpose was to demonstrate the UI/UX design with the Modern Botanical Atelier aesthetic, not to provide actual AI-powered predictions. The app should display mock toxicity analysis results without making any external API calls.

Key files:
- @services/geminiService.ts - Handles API calls to Gemini
- @App.tsx - Main application component that uses the service
- @vite.config.ts - Defines environment variable configuration
</context>

<requirements>
1. Remove all Gemini API integration code and dependencies
2. Replace with mock data that demonstrates the UI functionality
3. Mock data should show realistic toxicity analysis results including:
   - Toxicity classification (e.g., "High", "Moderate", "Low", "Uncertain")
   - Confidence scores (0-100%)
   - Explanations of the analysis
   - Recommendations
4. Maintain all existing UI/UX styling and components
5. Remove the @google/genai dependency from package.json
6. Remove API key configuration from vite.config.ts
7. Ensure the app works immediately without any environment variables
</requirements>

<implementation>
1. Modify `./services/geminiService.ts`:
   - Remove the Google GenAI import and initialization
   - Replace the analyzeToxicity function with a mock implementation
   - Generate deterministic mock results based on input parameters
   - Include varied, realistic responses for different compound types
   - Keep the same function signature and return type for compatibility

2. Update `./vite.config.ts`:
   - Remove the API_KEY and GEMINI_API_KEY environment variable definitions
   - Clean up any related configuration

3. Update `./package.json`:
   - Remove the @google/genai dependency

4. Verify `./App.tsx` requires no changes (should work with mock service)

Mock data strategy:
- Use compound name/type to generate deterministic but varied results
- Include realistic toxicity classifications and confidence levels
- Provide detailed explanations that reference the compound properties
- Make recommendations appropriate for each toxicity level
</implementation>

<output>
Modify the following files:
- `./services/geminiService.ts` - Replace with mock implementation
- `./vite.config.ts` - Remove API key configuration
- `./package.json` - Remove @google/genai dependency
</output>

<verification>
Before declaring complete:
1. Verify the app builds successfully: `npm run build`
2. Confirm no references to GEMINI_API_KEY remain in the codebase
3. Test that different compound inputs produce varied mock results
4. Ensure all UI components display correctly with mock data
5. Check that no API calls are made (no network requests to external services)
</verification>

<success_criteria>
- App runs without requiring any environment variables
- Mock toxicity analysis displays realistic, varied results
- All UI/UX functionality works identically to before
- No @google/genai dependency in package.json
- No API key configuration in vite.config.ts
- Build completes without errors
</success_criteria>
