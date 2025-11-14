# Security Policy

## Reporting a Vulnerability

We take security seriously. If you discover a security vulnerability in bmyCure4MM, please report it responsibly.

### How to Report

**DO NOT** create a public GitHub issue for security vulnerabilities.

Instead, please send an email to:
- **andreazedda@example.com**

Include in your report:
- Description of the vulnerability
- Steps to reproduce the issue
- Potential impact
- Suggested fix (if any)

### What to Expect

- **Acknowledgment**: We'll acknowledge receipt within 48 hours
- **Assessment**: We'll assess the vulnerability and determine severity
- **Fix**: We'll work on a fix and keep you updated on progress
- **Disclosure**: Once fixed, we'll coordinate public disclosure

### Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 1.0.x   | :white_check_mark: |
| < 1.0   | :x:                |

## Security Best Practices

### For Deployment

1. **Environment Variables**
   - Never commit `.env` files
   - Use strong, random `DJANGO_SECRET_KEY`
   - Set `DJANGO_DEBUG=0` in production

2. **Database**
   - Use PostgreSQL in production
   - Enable SSL connections
   - Regular backups

3. **HTTPS**
   - Always use HTTPS in production
   - Enable HSTS headers
   - Set secure cookie flags

4. **Dependencies**
   - Regularly update dependencies
   - Monitor security advisories
   - Use `pip-audit` or similar tools

5. **Access Control**
   - Implement proper authentication
   - Use role-based permissions
   - Enable two-factor authentication if possible

### For Development

1. **Local Setup**
   - Keep local `.env` files secure
   - Don't share credentials
   - Use test databases for development

2. **Code Review**
   - Review all security-related changes
   - Check for SQL injection vulnerabilities
   - Validate all user inputs

3. **Testing**
   - Write security tests
   - Test authentication flows
   - Verify permission checks

## Known Security Considerations

### Input Validation
- All forms use Django's built-in validation
- Additional validation in model `clean()` methods
- XSS prevention through template auto-escaping

### SQL Injection
- Protected by Django ORM
- No raw SQL queries without parameterization
- Queryset filtering uses safe methods

### CSRF Protection
- Enabled by default in Django
- All forms include CSRF tokens
- API endpoints use token authentication

### Authentication
- Built on Django's authentication system
- Passwords hashed with PBKDF2
- Session security configured

## Security Features

### Current Implementation
- ✅ CSRF protection
- ✅ XSS prevention
- ✅ SQL injection protection
- ✅ Secure password hashing
- ✅ Session security
- ✅ Environment-based secrets

### Recommended Additions
- [ ] Rate limiting
- [ ] Two-factor authentication
- [ ] Security headers (django-csp)
- [ ] IP whitelisting for admin
- [ ] Audit logging
- [ ] Intrusion detection

## Vulnerability Disclosure Timeline

1. **Day 0**: Vulnerability reported
2. **Day 1-2**: Acknowledgment sent
3. **Day 3-7**: Assessment and severity determination
4. **Day 7-30**: Fix development and testing
5. **Day 30**: Coordinated public disclosure

## Attribution

Security researchers who responsibly disclose vulnerabilities will be credited in:
- Release notes
- Security advisories
- Project documentation (with permission)

## Questions?

For security questions that aren't vulnerabilities, please:
- Open a GitHub Discussion
- Contact maintainers via email

Thank you for helping keep bmyCure4MM secure!
